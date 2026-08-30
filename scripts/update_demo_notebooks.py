from __future__ import annotations

import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
COMPUTATIONAL_COMMANDS = {
    "Figure2": ["--reproduce", "--step", "2"],
    "Figure3": ["--reproduce"],
    "Figure4": ["--reproduce"],
    "Figure5": ["--reproduce"],
    "FigureS3": ["--reproduce"],
}
EXTERNAL_DATA_FIGURES = {"FigureS3"}


FIGURE5_METRIC_NOTE = """## Boundary metrics

`dropdown_edge` is the mean signed raw ligand-count difference across the four signal-loss levels; it is not normalized by tool output size. `equal_edge` and `gradient_edge` are two-sample Kolmogorov-Smirnov D statistics:

`D = max_x |F_original(x) - F_perturbed(x)|`

The KS statistic ranges from 0 to 1. A larger value indicates a stronger distribution shift, but does not encode its direction.
"""


def markdown(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": [text]}


def code(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


def update_notebook(figure_dir: Path) -> None:
    config = json.loads((figure_dir / "figure_config.json").read_text(encoding="utf-8"))
    figure = config["figure"]
    notebook_path = figure_dir / f"{figure}_demo.ipynb"
    status = config["status"]
    notes = config["notes"]

    cells = [
        markdown(f"# {figure}\n\n{config['caption']}\n\n**Status:** `{status}`\n"),
        markdown(
            "## Result contract\n\n"
            f"{notes}\n\n"
            "`results/generated/` contains the latest script-generated SVG/CSV/JSON files. "
            "`results/published/` contains the submitted assembly and may include manual layout.\n"
        ),
    ]

    if figure == "Figure5":
        cells.append(markdown(FIGURE5_METRIC_NOTE))

    if figure in COMPUTATIONAL_COMMANDS:
        args = COMPUTATIONAL_COMMANDS[figure]
        args_literal = ", ".join(repr(item) for item in args)
        cells.extend(
            [
                markdown(
                    "## Reproduce\n\n"
                    "This cell calls the same figure shortcut used from the command line. "
                    "Figure2 runs its repository-contained CDF step; Figure3-5 use their canonical latest inputs. "
                    "FigureS3 reads precomputed h5ad files from `SPATIALCCCBENCH_DATA_ROOT`.\n"
                ),
                code(
                    "from pathlib import Path\n"
                    "import os\n"
                    "import subprocess\n"
                    "import sys\n\n"
                    f"figure = {figure!r}\n"
                    "script = Path.cwd() / 'run.py'\n"
                    "if not script.exists():\n"
                    "    script = Path('figures') / figure / 'run.py'\n"
                    f"command = [sys.executable, str(script), {args_literal}]\n"
                    + (
                        "data_root = Path(os.environ['SPATIALCCCBENCH_DATA_ROOT'])\n"
                        "command.extend(['--data-root', str(data_root)])\n"
                        if figure in EXTERNAL_DATA_FIGURES
                        else ""
                    )
                    +
                    "subprocess.run(command, check=True)\n"
                ),
                markdown(
                    "## Inspect\n\n"
                    "Review `results/result_manifest.json` for input scope and SHA-256 hashes of generated files.\n"
                ),
            ]
        )
    else:
        cells.extend(
            [
                markdown(
                    "## Inspect submitted result\n\n"
                    "This figure is reference/manual or has no complete standalone redraw. "
                    "The notebook intentionally does not claim to regenerate it.\n"
                ),
                code(
                    "from pathlib import Path\n\n"
                    f"figure = {figure!r}\n"
                    "figure_dir = Path.cwd()\n"
                    "if not (figure_dir / 'results').exists():\n"
                    "    figure_dir = Path('figures') / figure\n"
                    "figure_dir / 'results' / 'published' / f'{figure}.svg'\n"
                ),
            ]
        )

    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    notebook_path.write_text(json.dumps(notebook, indent=1) + "\n", encoding="utf-8")


def main() -> None:
    for figure_dir in sorted((REPO_ROOT / "figures").glob("Figure*")):
        if (figure_dir / "figure_config.json").exists():
            update_notebook(figure_dir)
            print(f"Updated {figure_dir.name}")


if __name__ == "__main__":
    main()
