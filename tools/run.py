#!/usr/bin/env python3
"""Install isolated environments and run SpatialCCC benchmark tools."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = Path(__file__).with_name("tool_registry.json")


class ToolError(RuntimeError):
    pass


def load_registry() -> dict[str, dict[str, Any]]:
    with REGISTRY_PATH.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != 1 or not isinstance(payload.get("tools"), dict):
        raise ToolError(f"Unsupported tool registry: {REGISTRY_PATH}")
    return payload["tools"]


def alias_map(tools: dict[str, dict[str, Any]]) -> dict[str, str]:
    aliases: dict[str, str] = {}
    for tool_id, spec in tools.items():
        aliases[tool_id.lower()] = tool_id
        for alias in spec.get("aliases", []):
            aliases[str(alias).lower()] = tool_id
    return aliases


def resolve_tool(name: str, tools: dict[str, dict[str, Any]]) -> tuple[str, dict[str, Any]]:
    canonical = alias_map(tools).get(name.lower())
    if canonical is None:
        choices = ", ".join(sorted(tools))
        raise ToolError(f"Unknown tool '{name}'. Available tools: {choices}")
    return canonical, tools[canonical]


def resolve_manager(requested: str, required: bool = True) -> str | None:
    candidates = [requested] if requested != "auto" else ["micromamba", "mamba", "conda"]
    for candidate in candidates:
        executable = shutil.which(candidate)
        if executable:
            return executable
        if candidate == "conda":
            discovered = discover_conda()
            if discovered:
                return discovered
    if required:
        raise ToolError("Conda, Mamba, or Micromamba was not found in PATH.")
    return None


def discover_conda() -> str | None:
    """Find Conda from CONDA_EXE or a Python installed below <conda>/envs/."""
    explicit = os.environ.get("CONDA_EXE")
    if explicit and Path(explicit).is_file():
        return explicit

    roots = [Path(sys.prefix), Path.home() / "miniconda3", Path.home() / "anaconda3"]
    executable_path = Path(sys.executable).resolve()
    roots.extend(executable_path.parents)
    seen: set[Path] = set()
    for root in roots:
        if root in seen:
            continue
        seen.add(root)
        candidates = (
            root / "Scripts" / "conda.exe",
            root / "bin" / "conda",
            root / "condabin" / "conda.bat",
        )
        for candidate in candidates:
            if candidate.is_file():
                return str(candidate)
    return None


def environment_names(manager: str) -> set[str]:
    completed = subprocess.run(
        [manager, "env", "list", "--json"],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    return {Path(path).name for path in payload.get("envs", [])}


def format_command(command: Iterable[str]) -> str:
    values = [str(value) for value in command]
    return subprocess.list2cmdline(values) if os.name == "nt" else shlex.join(values)


def manager_run_prefix(manager: str, env_name: str) -> list[str]:
    prefix = [manager, "run"]
    if Path(manager).stem.lower() == "conda":
        prefix.append("--no-capture-output")
    return [*prefix, "-n", env_name]


def run_process(command: list[str], *, dry_run: bool, log_path: Path | None = None) -> None:
    print(f"$ {format_command(command)}")
    if dry_run:
        return

    log_handle = None
    try:
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_handle = log_path.open("w", encoding="utf-8", newline="")
            log_handle.write(f"$ {format_command(command)}\n\n")

        process = subprocess.Popen(
            command,
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            if log_handle is not None:
                log_handle.write(line)
        return_code = process.wait()
        if return_code:
            raise ToolError(f"Command failed with exit code {return_code}.")
    finally:
        if log_handle is not None:
            log_handle.close()


def print_tools(tools: dict[str, dict[str, Any]]) -> None:
    headings = ("tool", "status", "environment", "expected output")
    rows = [
        (
            tool_id,
            spec["status"],
            spec.get("environment") or "-",
            spec.get("output") or "-",
        )
        for tool_id, spec in tools.items()
    ]
    widths = [max(len(headings[i]), *(len(str(row[i])) for row in rows)) for i in range(4)]
    print("  ".join(headings[i].ljust(widths[i]) for i in range(4)))
    print("  ".join("-" * width for width in widths))
    for row in rows:
        print("  ".join(str(row[i]).ljust(widths[i]) for i in range(4)))


def selected_tools(names: list[str], tools: dict[str, dict[str, Any]]) -> list[tuple[str, dict[str, Any]]]:
    if not names or names == ["all"]:
        return list(tools.items())
    selected: list[tuple[str, dict[str, Any]]] = []
    seen: set[str] = set()
    for name in names:
        tool_id, spec = resolve_tool(name, tools)
        if tool_id not in seen:
            selected.append((tool_id, spec))
            seen.add(tool_id)
    return selected


def install_tools(args: argparse.Namespace, tools: dict[str, dict[str, Any]]) -> None:
    resolved_manager = resolve_manager(args.manager, required=not args.dry_run)
    manager = resolved_manager or (
        "conda" if args.manager == "auto" else args.manager
    )
    existing = environment_names(manager) if resolved_manager else set()
    env_specs: dict[str, tuple[Path, list[str] | None]] = {}

    install_all = not args.tools or args.tools == ["all"]
    for tool_id, spec in selected_tools(args.tools, tools):
        if spec.get("platform_note"):
            print(f"[platform note] {tool_id}: {spec['platform_note']}")
        if install_all and spec["status"] == "unavailable":
            print(
                f"[skip] {tool_id}: dependency snapshot is not installed by default "
                "because its runner is unavailable"
            )
            continue
        env_name = spec.get("environment")
        env_file = spec.get("environment_file")
        if not env_name or not env_file:
            print(f"[skip] {tool_id}: {spec.get('reason', 'no managed environment')}")
            continue
        env_specs[env_name] = (REPO_ROOT / env_file, spec.get("post_install"))

    for env_name, (env_file, post_install) in env_specs.items():
        if not env_file.is_file():
            raise ToolError(f"Environment definition is missing: {env_file}")
        if env_name in existing and args.recreate:
            print(f"[remove] {env_name}")
            run_process(
                [manager, "env", "remove", "-n", env_name, "-y"],
                dry_run=args.dry_run,
            )
            existing.discard(env_name)
        if env_name in existing and not args.recreate:
            command = [manager, "env", "update", "-n", env_name, "-f", str(env_file), "--prune"]
            action = "update"
        else:
            command = [manager, "env", "create", "-f", str(env_file)]
            action = "create"
        print(f"[{action}] {env_name}")
        run_process(command, dry_run=args.dry_run)
        if post_install:
            script, *post_args = post_install
            post_command = [
                *manager_run_prefix(manager, env_name),
                "Rscript",
                str(REPO_ROOT / script),
                *post_args,
            ]
            print(f"[post-install] {env_name}")
            run_process(post_command, dry_run=args.dry_run)


def validate_input(path_value: str | None, label: str, dry_run: bool) -> Path:
    if not path_value:
        raise ToolError(f"{label} is required.")
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    path = path.resolve()
    if not dry_run and not path.is_file():
        raise ToolError(f"{label} does not exist: {path}")
    return path


def build_tool_command(
    tool_id: str,
    spec: dict[str, Any],
    args: argparse.Namespace,
    manager: str | None,
) -> tuple[list[str], Path]:
    if spec["status"] == "unavailable":
        raise ToolError(f"{tool_id} cannot be launched: {spec.get('reason', 'runner unavailable')}")

    runner = REPO_ROOT / spec["runner"]
    if not runner.is_file():
        raise ToolError(f"Runner is missing: {runner}")

    output_root = Path(args.output_root).expanduser()
    if not output_root.is_absolute():
        output_root = REPO_ROOT / output_root
    output_root = output_root.resolve()
    common = [
        "--analysis_dataset", args.dataset,
        "--LR_ref_path", str(validate_input(args.lr_db, "--lr-db", args.dry_run)),
        "--output_root", str(output_root),
    ]

    if spec["input_mode"] == "h5ad":
        common[0:0] = ["--adata_path", str(validate_input(args.adata, "--adata", args.dry_run))]
    elif spec["input_mode"] == "csv_pair":
        common[0:0] = [
            "--count_path", str(validate_input(args.count_csv, "--count-csv", args.dry_run)),
            "--meta_path", str(validate_input(args.meta_csv, "--meta-csv", args.dry_run)),
        ]

    runtime = spec["runtime"]
    executable = "python" if runtime == "python" else "Rscript"
    if args.active_environment:
        if runtime == "python":
            executable = sys.executable
        else:
            executable = shutil.which("Rscript") or "Rscript"
        command = [executable, str(runner)]
    else:
        assert manager is not None
        command = [*manager_run_prefix(manager, spec["environment"]), executable, str(runner)]

    command.extend(common)
    command.extend(str(value) for value in spec.get("default_args", []))
    command.extend(args.runner_arg or [])
    log_path = output_root / args.dataset / "_logs" / f"{tool_id}.log"
    return command, log_path


def run_tool(args: argparse.Namespace, tools: dict[str, dict[str, Any]]) -> None:
    tool_id, spec = resolve_tool(args.tool, tools)
    if spec.get("platform_note"):
        print(f"[platform note] {tool_id}: {spec['platform_note']}")
    if spec["status"] == "legacy":
        print(f"[warning] {tool_id} uses the manuscript-era R API; Linux is recommended.")
    manager = None if args.active_environment else resolve_manager(args.manager, required=not args.dry_run)

    if manager and not args.dry_run and not args.skip_environment_check:
        installed = environment_names(manager)
        if spec.get("environment") not in installed:
            raise ToolError(
                f"Environment '{spec.get('environment')}' is not installed. "
                f"Run: python tools/run.py install {tool_id}"
            )

    fallback_manager = "conda" if args.manager == "auto" else args.manager
    command, log_path = build_tool_command(tool_id, spec, args, manager or fallback_manager)
    run_process(command, dry_run=args.dry_run, log_path=log_path)
    if not args.dry_run:
        output_root = Path(args.output_root).expanduser()
        if not output_root.is_absolute():
            output_root = REPO_ROOT / output_root
        expected = output_root.resolve() / args.dataset / spec["output"]
        if not expected.exists():
            raise ToolError(f"Runner finished, but expected output was not created: {expected}")
        print(f"Output: {expected}")
        print(f"Log: {log_path}")


def reproduce_tool(args: argparse.Namespace, tools: dict[str, dict[str, Any]]) -> None:
    if args.active_environment:
        print("[skip install] --active-environment was requested.")
    else:
        install_args = argparse.Namespace(
            tools=[args.tool],
            manager=args.manager,
            recreate=args.recreate,
            dry_run=args.dry_run,
        )
        install_tools(install_args, tools)
    run_tool(args, tools)


def doctor(args: argparse.Namespace, tools: dict[str, dict[str, Any]]) -> None:
    manager = resolve_manager(args.manager, required=False)
    installed: set[str] = set()
    if manager:
        try:
            installed = environment_names(manager)
        except (subprocess.SubprocessError, json.JSONDecodeError) as exc:
            print(f"[fail] environment manager: {exc}")
    print(f"Repository: {REPO_ROOT}")
    print(f"Environment manager: {manager or 'not found'}")

    failures = 0
    for tool_id, spec in selected_tools(args.tools, tools):
        runner = REPO_ROOT / spec["runner"] if spec.get("runner") else None
        env_file = REPO_ROOT / spec["environment_file"] if spec.get("environment_file") else None
        runner_ok = bool(runner and runner.is_file())
        env_file_ok = bool(env_file and env_file.is_file())
        env_ok = bool(spec.get("environment") and spec["environment"] in installed)
        if spec["status"] == "unavailable":
            print(
                f"[unavailable runner] {tool_id}: "
                f"env-file={'yes' if env_file_ok else 'no'}, "
                f"installed={'yes' if env_ok else 'no'}; {spec.get('reason', '')}"
            )
            failures += int(not env_file_ok)
            continue
        print(
            f"[{'ok' if runner_ok and env_file_ok else 'fail'}] {tool_id}: "
            f"runner={'yes' if runner_ok else 'no'}, env-file={'yes' if env_file_ok else 'no'}, "
            f"installed={'yes' if env_ok else 'no'}"
        )
        failures += int(not runner_ok or not env_file_ok)
    if failures:
        raise ToolError(f"Doctor found {failures} incomplete tool definition(s).")


def run_batch(args: argparse.Namespace, tools: dict[str, dict[str, Any]]) -> None:
    batch_path = validate_input(args.jobs, "--jobs", args.dry_run)
    with batch_path.open(encoding="utf-8-sig", newline="") as handle:
        jobs = list(csv.DictReader(handle))
    if not jobs:
        raise ToolError(f"No jobs found in {batch_path}")

    for index, row in enumerate(jobs, start=1):
        if not row.get("tool") or not row.get("dataset"):
            raise ToolError(f"Batch row {index} must define tool and dataset.")
        print(f"\n[{index}/{len(jobs)}] {row['tool']} on {row['dataset']}")
        run_args = argparse.Namespace(
            tool=row["tool"],
            dataset=row["dataset"],
            adata=row.get("adata") or None,
            lr_db=row.get("lr_db") or None,
            count_csv=row.get("count_csv") or None,
            meta_csv=row.get("meta_csv") or None,
            output_root=row.get("output_root") or args.output_root,
            manager=args.manager,
            active_environment=args.active_environment,
            skip_environment_check=args.skip_environment_check,
            runner_arg=[],
            dry_run=args.dry_run,
        )
        run_tool(run_args, tools)


def add_manager_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manager", choices=("auto", "conda", "mamba", "micromamba"), default="auto")


def add_run_arguments(parser: argparse.ArgumentParser, include_tool: bool = True) -> None:
    if include_tool:
        parser.add_argument("tool")
    parser.add_argument("--dataset", required=True, help="Dataset label used below the output root.")
    parser.add_argument("--adata", help="Input h5ad for Python tools.")
    parser.add_argument("--lr-db", required=True, help="CellChat-style ligand-receptor CSV.")
    parser.add_argument("--count-csv", help="Spot-by-gene count CSV for R tools.")
    parser.add_argument("--meta-csv", help="Spot metadata CSV for R tools.")
    parser.add_argument("--output-root", default=str(REPO_ROOT / "result"))
    parser.add_argument("--active-environment", action="store_true", help="Run in the active environment instead of conda run.")
    parser.add_argument("--skip-environment-check", action="store_true")
    parser.add_argument("--runner-arg", action="append", default=[], help="Additional runner argument; repeat for each token.")
    parser.add_argument("--dry-run", action="store_true")
    add_manager_argument(parser)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("list", help="List registered benchmark tools.")

    doctor_parser = subparsers.add_parser("doctor", help="Check runners, environment files, and installed environments.")
    doctor_parser.add_argument("tools", nargs="*", default=["all"])
    add_manager_argument(doctor_parser)

    install_parser = subparsers.add_parser("install", help="Create or update isolated tool environments.")
    install_parser.add_argument("tools", nargs="+", help="Tool names or 'all'.")
    install_parser.add_argument("--recreate", action="store_true", help="Create rather than update existing environments.")
    install_parser.add_argument("--dry-run", action="store_true")
    add_manager_argument(install_parser)

    run_parser = subparsers.add_parser("run", help="Run one benchmark tool.")
    add_run_arguments(run_parser)

    reproduce_parser = subparsers.add_parser("reproduce", help="Install one environment and then run its tool.")
    add_run_arguments(reproduce_parser)
    reproduce_parser.add_argument("--recreate", action="store_true")

    batch_parser = subparsers.add_parser("batch", help="Run jobs from a CSV manifest.")
    batch_parser.add_argument("--jobs", required=True)
    batch_parser.add_argument("--output-root", default=str(REPO_ROOT / "result"))
    batch_parser.add_argument("--active-environment", action="store_true")
    batch_parser.add_argument("--skip-environment-check", action="store_true")
    batch_parser.add_argument("--dry-run", action="store_true")
    add_manager_argument(batch_parser)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    tools = load_registry()
    try:
        if args.command == "list":
            print_tools(tools)
        elif args.command == "doctor":
            doctor(args, tools)
        elif args.command == "install":
            install_tools(args, tools)
        elif args.command == "run":
            run_tool(args, tools)
        elif args.command == "reproduce":
            reproduce_tool(args, tools)
        elif args.command == "batch":
            run_batch(args, tools)
        else:
            parser.error(f"Unsupported command: {args.command}")
    except (ToolError, subprocess.CalledProcessError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
