import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import math
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from matplotlib.patches import Rectangle
from matplotlib import patches
from matplotlib.path import Path as MplPath


NOISE_ORDER = ["dropout", "non_specific", "up", "lack", "overlap", "offset"]
NOISE_LABELS = {
    "dropout": "Random dropout",
    "non_specific": "Non-specific transcription noise",
    "up": "Random interpolation",
    "lack": "Tissue section lack",
    "overlap": "Tissue overlap",
    "offset": "Tissue fracture",
}
NOISE_METRIC_STYLE = {
    "precision": {"label": "Precision", "color": "#2E86C1", "marker": "o"},
    "recall": {"label": "Recall", "color": "#E74C3C", "marker": "s"},
    "F1": {"label": "F1 score", "color": "#F39C12", "marker": "^"},
}


def _polar_to_axes(theta, radius, max_radius=1.0):
    scaled = 0.5 * radius / max_radius
    return 0.5 + scaled * math.cos(theta), 0.5 + scaled * math.sin(theta)


def _add_noise_outer_segment(ax, theta_center, label):
    gap = math.radians(8)
    span = 2 * math.pi / len(NOISE_ORDER) - gap
    theta_start = theta_center - span / 2
    theta_end = theta_center + span / 2
    outer_thetas = np.linspace(theta_start, theta_end, 40)
    inner_thetas = np.linspace(theta_end, theta_start, 40)
    vertices = [_polar_to_axes(float(theta), 1.18) for theta in outer_thetas]
    vertices += [_polar_to_axes(float(theta), 1.05) for theta in inner_thetas]
    vertices.append(vertices[0])
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(vertices) - 2) + [MplPath.CLOSEPOLY]
    ax.add_patch(
        patches.PathPatch(
            MplPath(vertices, codes),
            transform=ax.transAxes,
            facecolor="#363636",
            edgecolor="none",
            clip_on=False,
            zorder=30,
        )
    )

    cleaned = " ".join(label.split())
    char_count = len(cleaned)
    text_span = min(span * 0.90, math.radians(max(18, char_count * 3.0)))
    upper_half = math.sin(theta_center) >= 0
    if char_count == 1:
        text_thetas = [theta_center]
    elif upper_half:
        text_thetas = np.linspace(theta_center + text_span / 2, theta_center - text_span / 2, char_count)
    else:
        text_thetas = np.linspace(theta_center - text_span / 2, theta_center + text_span / 2, char_count)

    for char, theta in zip(cleaned, text_thetas):
        theta = float(theta)
        x, y = _polar_to_axes(theta, 1.115)
        rotation = math.degrees(theta) - 90 if upper_half else math.degrees(theta) + 90
        ax.text(
            x,
            y,
            char,
            transform=ax.transAxes,
            ha="center",
            va="center",
            rotation=rotation,
            rotation_mode="anchor",
            fontsize=4.2,
            color="white",
            clip_on=False,
            zorder=31,
        )


def prepare_noise_metrics(df, tool, noise_order=NOISE_ORDER):
    """Validate and order one tool's processed noise metrics."""
    required = {"precision", "recall", "F1", "tool", "noise"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Noise table is missing columns: {sorted(missing)}")

    selected = df.loc[df["tool"].astype(str) == str(tool)].copy()
    if selected["noise"].duplicated().any():
        duplicated = selected.loc[selected["noise"].duplicated(), "noise"].tolist()
        raise ValueError(f"Duplicate noise rows for {tool}: {duplicated}")
    selected = selected.set_index("noise").reindex(noise_order)
    if selected[["precision", "recall", "F1"]].isna().any().any():
        missing_noise = selected.index[selected["F1"].isna()].tolist()
        raise ValueError(f"Missing noise metrics for {tool}: {missing_noise}")
    values = selected[["precision", "recall", "F1"]].astype(float)
    if ((values < 0) | (values > 1)).any().any():
        raise ValueError(f"Noise metrics for {tool} must be between 0 and 1")
    return selected


def plot_noise_resistance_spider(df, tool, save_path, noise_order=NOISE_ORDER):
    """Render one reference-style noise-resistance radar plot as SVG."""
    selected = prepare_noise_metrics(df, tool, noise_order=noise_order)
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "svg.fonttype": "none",
        }
    )
    angles = np.linspace(0, 2 * np.pi, len(noise_order), endpoint=False)
    closed_angles = np.concatenate([angles, angles[:1]])

    fig, ax = plt.subplots(figsize=(3.35, 3.45), subplot_kw={"polar": True})
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    fig.subplots_adjust(top=0.80, bottom=0.08, left=0.10, right=0.90)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_xticks(angles)
    ax.set_xticklabels([""] * len(angles))
    ax.set_yticks([0.2, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.6", "0.8", "1"], fontsize=7, color="black")
    ax.set_rlabel_position(30)
    ax.grid(color="#333333", linestyle="--", linewidth=0.45, alpha=0.75)
    ax.spines["polar"].set_visible(False)

    for metric, style in NOISE_METRIC_STYLE.items():
        values = selected.loc[noise_order, metric].to_numpy(dtype=float)
        closed_values = np.concatenate([values, values[:1]])
        ax.plot(
            closed_angles,
            closed_values,
            color=style["color"],
            linewidth=1.55,
            marker=style["marker"],
            markersize=4.2,
            markerfacecolor=style["color"],
            markeredgecolor="white",
            markeredgewidth=0.35,
            zorder=20,
        )
        ax.fill(closed_angles, closed_values, color=style["color"], alpha=0.16, zorder=10)

    display_angles = np.pi / 2 - angles
    for theta, noise in zip(display_angles, noise_order):
        _add_noise_outer_segment(ax, float(theta), NOISE_LABELS[noise])

    title_size = 11 if len(str(tool)) > 20 else 13
    ax.set_title(str(tool), fontsize=title_size, fontweight="bold", pad=28)
    output_format = save_path.suffix.lower().lstrip(".") or "svg"
    fig.savefig(save_path, format=output_format, dpi=300, bbox_inches="tight", pad_inches=0.14, facecolor="white")
    plt.close(fig)
    return save_path


def plot_noise_metric_legend(save_path, ncol=3):
    """Save a standalone Precision/Recall/F1 score legend as SVG."""
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    handles = [
        Line2D(
            [0],
            [0],
            color=style["color"],
            marker=style["marker"],
            linewidth=1.8,
            markersize=5.5,
            markerfacecolor=style["color"],
            markeredgecolor="white",
            markeredgewidth=0.4,
            label=style["label"],
        )
        for style in NOISE_METRIC_STYLE.values()
    ]
    fig, ax = plt.subplots(figsize=(4.5, 0.55))
    ax.axis("off")
    ax.legend(handles=handles, loc="center", ncol=ncol, frameon=False, fontsize=9, handlelength=2.4)
    output_format = save_path.suffix.lower().lstrip(".") or "svg"
    fig.savefig(save_path, format=output_format, dpi=300, bbox_inches="tight", pad_inches=0.03, transparent=True)
    plt.close(fig)
    return save_path


def plot_all_noise_resistance(df, out_dir):
    """Render one SVG per tool plus a standalone metric legend."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for tool in df["tool"].drop_duplicates().astype(str):
        paths.append(plot_noise_resistance_spider(df, tool, out_dir / f"{tool}_spider_plot.svg"))
    legend_path = plot_noise_metric_legend(out_dir / "precision_recall_f1_legend.svg")
    return paths, legend_path

def cal_normal_result(metrics_df):
    metrics_df['Precision_norm']= metrics_df.groupby('Tool')['Precision'].transform(lambda x: (x - x.min()) / (x.max() - x.min()))
    metrics_df['Recall_norm']= metrics_df.groupby('Tool')['Recall'].transform(lambda x: (x - x.min()) / (x.max() - x.min()))
    metrics_df['F1 Score_norm']= metrics_df.groupby('Tool')['F1 Score'].transform(lambda x: (x - x.min()) / (x.max() - x.min()))
    metrics_df['Precision_norm']=metrics_df['Precision_norm']*100
    metrics_df['Precision_mean']= metrics_df.groupby('Tool')['Precision_norm'].transform(lambda x: (x.sum()/len(x)))
    metrics_df['Recall_norm']=metrics_df['Recall_norm']*100
    metrics_df['Recall_mean']= metrics_df.groupby('Tool')['Recall_norm'].transform(lambda x: (x.sum()/len(x)))
    metrics_df['F1 Score_norm']=metrics_df['F1 Score_norm']*100
    metrics_df['F1 Score_mean']= metrics_df.groupby('Tool')['F1 Score_norm'].transform(lambda x: (x.sum()/len(x)))
    #metrics_df['Recall_mean']= metrics_df.groupby('Tool')['Recall'].transform(lambda x: (x - x.min()) / (x.max() - x.min()))
    #metrics_df['F1 Score_mean']= metrics_df.groupby('Tool')['F1 Score'].transform(lambda x: (x - x.min()) / (x.max() - x.min()))
    df = metrics_df

    return df

def plot_network(tool, df, cell_color_mapping, node_order=None,analysis_target="Unknown", base_color="#3498db"):
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['font.weight'] = 'normal'  # 设置为正常不加粗


    if 'label' not in df.columns:
        raise ValueError("DataFrame must contain 'label' column")


    all_nodes = set()
    for idx, row in df.iterrows():
        source, target = extract_cell_types(row['Cell_Interaction'])
        if source:
            all_nodes.add(source)
        if target:
            all_nodes.add(target)

    G = nx.DiGraph()
    G.add_nodes_from(all_nodes)


    edges = []
    edge_labels = {}
    for idx, row in df.iterrows():
        source, target = extract_cell_types(row['Cell_Interaction'])
        if source and target:
            edges.append((source, target, row['Ligand_Receptor']))

            edge_key = (source, target)

            if edge_key not in edge_labels:
                edge_labels[edge_key] = row['label']
            else:

                edge_labels[edge_key] = min(edge_labels[edge_key], row['label'])


    edge_weights = {}
    for source, target, lr_pair in edges:
        key = (source, target)
        if key in edge_weights:
            edge_weights[key] += 1
        else:
            edge_weights[key] = 1


    for (source, target), weight in edge_weights.items():

        G.add_edge(source, target, weight=weight, label=edge_labels.get((source, target), 0))

    if node_order:

        missing_nodes = set(node_order) - set(all_nodes)
        if missing_nodes:
            print(f"Added missing nodes: {missing_nodes}")

            G.add_nodes_from(missing_nodes)

            all_nodes = set(G.nodes())
    else:

        node_order = sorted(all_nodes)

    plt.figure(figsize=(14, 14))

    pos = nx.circular_layout(G, scale=1.0, center=(0, 0), dim=2)

    if node_order:
        ordered_pos = {}
        angle_step = 2 * np.pi / len(node_order)

        for i, node in enumerate(node_order):
            if node in pos:
                x, y = pos[node]
                angle = np.arctan2(y, x)
                new_angle = i * angle_step
                ordered_pos[node] = (np.cos(new_angle), np.sin(new_angle))
            else:
                ordered_pos[node] = (np.cos(i * angle_step), np.sin(i * angle_step))

        pos = ordered_pos

    gradient_light = generate_gradient_simple(base_color, 5, 'to light')

    weights_label0 = []
    weights_label1 = []

    for u, v in G.edges():
        edge_label = G[u][v].get('label', 0)
        weight = G[u][v]['weight']

        if edge_label == 1:
            weights_label1.append(weight)
        else:
            weights_label0.append(weight)

    min_weight0 = min(weights_label0) if weights_label0 else 1
    max_weight0 = max(weights_label0) if weights_label0 else 1
    min_weight1 = min(weights_label1) if weights_label1 else 1
    max_weight1 = max(weights_label1) if weights_label1 else 1

    nodes = list(G.nodes())

    node_colors = []
    for node in nodes:
        if node in cell_color_mapping:
            node_colors.append(cell_color_mapping[node])
        else:
            import hashlib
            hash_obj = hashlib.md5(node.encode())
            hash_int = int(hash_obj.hexdigest()[:8], 16)
            hue = (hash_int % 360) / 360.0
            import colorsys
            rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
            hex_color = f"#{int(rgb[0]*255):02x}{int(rgb[1]*255):02x}{int(rgb[2]*255):02x}"
            node_colors.append(hex_color)
            print(f"Warning: node '{node}' isnt in color list added: {hex_color}")

    edge_colors = []
    for u, v in G.edges():
        weight = G[u][v]['weight']
        if weights_label0:
            if weight > np.percentile(weights_label0, 80):
                edge_colors.append(gradient_light[4])
            elif weight > np.percentile(weights_label0, 60):
                edge_colors.append(gradient_light[3])
            elif weight > np.percentile(weights_label0, 40):
                edge_colors.append(gradient_light[2])
            elif weight > np.percentile(weights_label0, 20):
                edge_colors.append(gradient_light[1])
            else:
                edge_colors.append(gradient_light[0])
        else:
            edge_colors.append(base_color)

    node_degrees = dict(G.degree())
    max_degree = max(node_degrees.values()) if node_degrees else 1

    nx.draw_networkx_nodes(G, pos,
                           node_size=500,
                           node_color=node_colors,
                           alpha=0.8,
                           edgecolors=base_color,
                           linewidths=1.5)

    nx.draw_networkx_labels(G, pos,
                            font_size=9,
                            font_weight='bold',
                            font_color='black')


    for i, (u, v) in enumerate(G.edges()):
        edge_label = G[u][v].get('label', 0)
        weight = G[u][v]['weight']

        if edge_label == 1:
            width = 0.5 + 3 * (weight / max_weight1)
        else:
            width = 0.5 + 3 * (weight / max_weight0)

        nx.draw_networkx_edges(G, pos,
                              edgelist=[(u, v)],
                              width=width,
                              edge_color=edge_colors[i],
                              alpha=0.7,
                              arrows=True,
                              arrowsize=20,
                              arrowstyle='->',
                              connectionstyle='arc3,rad=0.1')

    legend_elements = []
    used_cell_types = set(nodes)

    for cell_type in used_cell_types:
        if cell_type in cell_color_mapping:
            color = cell_color_mapping[cell_type]
            legend_elements.append(plt.Line2D([0], [0],
                                            marker='o',
                                            color='w',
                                            markerfacecolor=color,
                                            markersize=8,
                                            label=cell_type))

    legend_elements.append(plt.Line2D([0], [0],
                                      color='#00FF00',
                                      lw=2,
                                      label='Label = 1'))

    legend_elements.append(plt.Line2D([0], [0],
                                      color=gradient_light[2],
                                      lw=2,
                                      label='Label = 0'))

    plt.text(0.5, 0.98, f'CCC Network of {tool}', fontweight='bold',
         transform=plt.gca().transAxes, fontsize=20,
         horizontalalignment='center',
         verticalalignment='top',
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))

    plt.axis('off')
    plt.tight_layout()


    plt.savefig(f'./figure/CCC net_work({tool.replace("*","without sp")}_{analysis_target}).svg', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')

    plt.show()

def plot_bar_scatter(df, parameter, data_name="",show_legend=True,figsize=(8, 8)):

    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['font.weight'] = 'normal'  # 设置为正常不加粗

    fig, ax = plt.subplots(figsize=figsize, dpi=120)

    fig.patch.set_facecolor('#FFFFFF')
    ax.set_facecolor('#FFFFFF')

    rectangle_width = 0.7
    tool_stats = df.groupby('Tool')[parameter].agg(['mean', 'std', 'count']).sort_values('mean', ascending=False)
    tools = tool_stats.index
    means = tool_stats['mean'].values
    errors = tool_stats['std'].values

    tool_colors = {
        'Squidpy': '#1F77B4',
        'SpatialDM': '#17BECF',
        'Baseline_1': '#006400',
        'stLearn': '#D62728',
        'Giotto': '#7F7F7F',
        'Baseline_2': '#8C564B',
        'stLearn*': '#E377C2',
        'CellAgentChat': '#9467BD',
        'SpaTalk': '#BCBD22',
        'COMMOT': '#FF7F0E',
    }

    for i, tool in enumerate(tools):
        tool_data = df[df['Tool'] == tool]
        tool_data = tool_data[parameter]
        q1s = np.percentile(tool_data, 25)
        q3s = np.percentile(tool_data, 75)
        iqr=q3s-q1s
        mean_val = np.mean(tool_data)
        std_val = np.std(tool_data)
        if tool in tools:
            tool_color=tool_colors[tool]
        else:
            tool_color="#FF00FF"
        rect = Rectangle(xy=(i - rectangle_width/2, 0),
                         width=rectangle_width,
                         height=mean_val,
                         facecolor=tool_color,
                         edgecolor=tool_color,
                         alpha=0.85,
                         linewidth=1.5,
                         zorder=3)
        ax.add_patch(rect)

        ax.plot([i, i], [mean_val-std_val, mean_val+std_val],
        color=tool_colors[tool], linewidth=2.5, alpha=0.85, zorder=3)

    for i, tool in enumerate(tools):
        tool_data = df[df['Tool'] == tool]
        tool_data = tool_data[parameter]
        color = tool_colors[tool]

        x_pos = np.random.normal(i, 0.15, size=len(tool_data))
        ax.scatter(x_pos, tool_data, color=color, alpha=0.65, s=10,
                   edgecolor='w', linewidth=0.5, zorder=4,
                   label='Individual cell_pairs' if i == 0 else "")


    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels(tools, rotation=45, ha='right', rotation_mode='anchor')  # 设置刻度标签为工具名

    ax.set_ylabel(f"{parameter}", fontsize=12, fontweight='bold', labelpad=10)
    ax.set_xlabel('Tool', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title(f'{parameter}\n({data_name})', fontsize=14, pad=20, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#CCCCCC')
    ax.spines['bottom'].set_color('#CCCCCC')

    ax.grid(axis='y', linestyle='--', alpha=0.7, color='#FFFFFF', zorder=1)

    y_min = 0
    y_max = df[parameter].max() +0.1
    ax.set_ylim(y_min, y_max)

    if show_legend:
        legend_elements = [
            Patch(facecolor='#CCCCCC', edgecolor='w', label=f'Average {parameter}'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#CCCCCC',
                   markersize=8, label='Individual cell_pairs', markeredgecolor='w')
        ]

        for tool, color in tool_colors.items():
            if tool in tools:
                legend_elements.append(
                    Patch(facecolor=color, edgecolor='w', label=tool)
                )

        legend = ax.legend(handles=legend_elements, loc='upper left',
                           bbox_to_anchor=(1.02, 1), borderaxespad=0.,
                           frameon=True, fancybox=True,
                           framealpha=0.9, edgecolor='#DDDDDD')
        legend.get_frame().set_facecolor('#FFFFFF')

    plt.tight_layout()
    plt.subplots_adjust(top=0.92, bottom=0.15, right=0.85 if show_legend else 0.95)

    plt.savefig(f'./figure/Precision/({data_name}_{parameter}).svg', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')

    plt.show()

def extract_cell_types(interaction_str):

    if '-' in interaction_str:
        parts = interaction_str.split('-')
        source = parts[0]
        target = parts[-1]
        return source, target
    else:
        return None, None

def generate_gradient_simple(base_hex, num_colors=5, direction='to dark', ratio=0.3):

    hex_color = base_hex.lstrip('#')
    base_r = int(hex_color[0:2], 16)
    base_g = int(hex_color[2:4], 16)
    base_b = int(hex_color[4:6], 16)

    gradient_colors = []

    for i in range(num_colors): # 'to dark'
        if direction == 'to dark':

            factor = 1 - (i / (num_colors - 1)) * ratio
            new_r = int(base_r * factor)
            new_g = int(base_g * factor)
            new_b = int(base_b * factor)
        else:  # 'to light'

            factor = (i / (num_colors - 1)) * ratio
            new_r = int(base_r + (255 - base_r) * factor)
            new_g = int(base_g + (255 - base_g) * factor)
            new_b = int(base_b + (255 - base_b) * factor)

        new_r = max(0, min(255, new_r))
        new_g = max(0, min(255, new_g))
        new_b = max(0, min(255, new_b))

        hex_color = f"#{new_r:02x}{new_g:02x}{new_b:02x}"
        gradient_colors.append(hex_color)

    return gradient_colors

def plot_auto_corr_bar(df_percentage,gini_index=None,colors=["pink", 'lightblue', 'red', 'lightyellow', '#90EE90'],title="Per Genes Mode",M="1",show_legend=True,y_label='Analysis Tools'):

    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['font.weight'] = 'normal'  # 设置为正常不加粗

    if 'sum' in df_percentage.columns:
        df_percentage = df_percentage.drop('sum', axis=1)

    plt.figure(figsize=(5, 10))
    ax = plt.subplot(111)

    bars = df_percentage.plot(kind='bar', stacked=True, ax=ax, color=colors,
                               linewidth=0.5, width=0.85)

    ax_title=f'{title}\n (Gini_index>{gini_index})'
    ax.set_title(ax_title,
                 fontsize=16, fontweight='bold', pad=20, color='#2c3e50')
    ax.set_xlabel(y_label, fontsize=13, fontweight='semibold', labelpad=10, color='#34495e')
    ax.set_ylabel('Percentage Distribution', fontsize=13, fontweight='semibold', labelpad=10, color='#34495e')

    ax.tick_params(axis='x', rotation=75, labelsize=11, colors='#2c3e50')
    ax.tick_params(axis='y', labelsize=11, colors='#2c3e50')

    ax.grid(True, axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)


    ax.set_ylim(0, 1)


    if show_legend:

        legend = ax.legend(title='Interaction Categories', bbox_to_anchor=(1.05, 1),
                           loc='upper left', frameon=True, fancybox=True,
                           shadow=True, framealpha=0.9, title_fontsize=12, fontsize=10)
        legend.get_frame().set_facecolor('#f8f9fa')
        legend.get_frame().set_edgecolor('#dee2e6')
    else:

        ax.legend().set_visible(False)


    for container in ax.containers:

        ax.bar_label(container, labels=[f'{v:.1f}%' if v > 5 else '' for v in container.datavalues],
                     label_type='center', fontsize=9, color='white', fontweight='bold')


    for spine in ax.spines.values():
        spine.set_color('#bdc3c7')
        spine.set_linewidth(0.8)

    ax.set_facecolor('#fafafa')
    plt.gcf().patch.set_facecolor('white')


    if show_legend:

        plt.tight_layout()
        plt.subplots_adjust(right=0.85)
    else:

        plt.tight_layout()

    plt.savefig(f'./figure/Auto_correlation/{title}_{M}.svg', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
    plt.show()


def get_intersection_color(intersection_tuple, mapping, tool_colors):
    included_tools = []
    for i, included in enumerate(intersection_tuple):
        if included:
            for tool, idx in mapping.items():
                if idx == i:
                    included_tools.append(tool)
                    break

    if not included_tools:
        return '#808080'

    first_tool = included_tools[0]
    return tool_colors.get(first_tool, '#808080')

def plot_kde_cdf_comparison(df, variable_names, variable_labels,
                           colors=None, line_styles=None, line_widths=None,
                           alpha_values=None, save_path='kde_smooth_cdf.svg',
                           title='Cumulative Distribution Functions (KDE Smoothed)',
                           xlabel='Value', ylabel='Cumulative Probability',
                           figsize=(12, 8), dpi=300, transparent=True,
                           n_smooth_points=2000, kde_bw='scott',
                           show_data_points=False, point_size=20, point_alpha=0.2,
                           show_legend=True, legend_loc='lower right',
                           grid_alpha=0.3, xlim=None, ylim=(0, 1.05)):
    """
    使用核密度估计(KDE)绘制平滑的CDF对比图

    参数:
    ----------
    df : pandas DataFrame
        包含数据的数据框
    variable_names : list
        要绘制的变量名列表
    variable_labels : list
        对应变量的标签列表（用于图例）
    colors : list, optional
        线条颜色列表
    line_styles : list, optional
        线条样式列表
    line_widths : list, optional
        线条宽度列表
    alpha_values : list, optional
        透明度列表
    save_path : str, optional
        保存图片的路径
    title : str, optional
        图表标题
    xlabel : str, optional
        x轴标签
    ylabel : str, optional
        y轴标签
    figsize : tuple, optional
        图表大小
    dpi : int, optional
        图片分辨率
    transparent : bool, optional
        是否保存为透明背景
    n_smooth_points : int, optional
        平滑曲线的点数
    kde_bw : str or float, optional
        KDE带宽方法
    show_data_points : bool, optional
        是否显示原始数据点
    point_size : int, optional
        数据点大小
    point_alpha : float, optional
        数据点透明度
    show_legend : bool, optional
        是否显示图例
    legend_loc : str, optional
        图例位置
    grid_alpha : float, optional
        网格透明度
    xlim : tuple, optional
        x轴范围
    ylim : tuple, optional
        y轴范围

    返回:
    ----------
    fig : matplotlib.figure.Figure
        图表对象
    ax : matplotlib.axes.Axes
        坐标轴对象
    """

    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)

    # 设置颜色
    if colors is None:
        colors = DEFAULT_COLORS[:len(variable_names)]
    elif len(colors) < len(variable_names):
        colors = [colors[i % len(colors)] for i in range(len(variable_names))]

    # 设置线条样式
    if line_styles is None:
        line_styles = ['-'] * len(variable_names)
    elif len(line_styles) < len(variable_names):
        line_styles = [line_styles[i % len(line_styles)] for i in range(len(variable_names))]

    # 设置线条宽度
    if line_widths is None:
        line_widths = [3] * len(variable_names)
    elif isinstance(line_widths, (int, float)):
        line_widths = [line_widths] * len(variable_names)
    elif len(line_widths) < len(variable_names):
        line_widths = [line_widths[i % len(line_widths)] for i in range(len(variable_names))]

    # 设置透明度
    if alpha_values is None:
        alpha_values = [0.9] * len(variable_names)
    elif isinstance(alpha_values, (int, float)):
        alpha_values = [alpha_values] * len(variable_names)
    elif len(alpha_values) < len(variable_names):
        alpha_values = [alpha_values[i % len(alpha_values)] for i in range(len(variable_names))]

    # 存储所有曲线的x轴范围
    all_x_min = []
    all_x_max = []

    # 绘制每条KDE平滑的CDF曲线
    for idx, (col_name, label) in enumerate(zip(variable_names, variable_labels)):
        values = df[col_name]

        # 获取原始经验CDF点
        sorted_values = np.sort(values)
        y_empirical = np.arange(1, len(sorted_values) + 1) / len(sorted_values)

        # 绘制原始数据点（可选）
        if show_data_points and len(sorted_values) < 500:  # 数据点太多时不显示
            ax.scatter(sorted_values, y_empirical, s=point_size,
                      alpha=point_alpha, color=colors[idx],
                      edgecolors='white', linewidths=0.5, zorder=5,
                      label=f'{label} (Data Points)')

        # 使用KDE平滑CDF
        xs_smooth, ys_smooth = kde_smooth_cdf(
            values, n_points=n_smooth_points, bw_method=kde_bw
        )

        # 存储x轴范围
        all_x_min.append(xs_smooth.min())
        all_x_max.append(xs_smooth.max())

        # 确保CDF在[0,1]范围内
        ys_smooth = np.clip(ys_smooth, 0, 1)

        # 绘制KDE平滑曲线
        ax.plot(xs_smooth, ys_smooth,
               linewidth=line_widths[idx],
               linestyle=line_styles[idx],
               label=label,
               color=colors[idx],
               alpha=alpha_values[idx],
               zorder=10)

    # 设置图表属性
    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    if show_legend:
        ax.legend(fontsize=11, frameon=False, loc=legend_loc)

    # 美化坐标轴
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

    # 设置坐标轴范围
    if xlim is None:
        ax.set_xlim(min(all_x_min), max(all_x_max))
    else:
        ax.set_xlim(xlim)

    ax.set_ylim(ylim)

    # 添加网格
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=grid_alpha)
    ax.set_axisbelow(True)

    plt.tight_layout()

    # 保存图片
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=transparent)
    print(f"KDE平滑CDF曲线图已保存为: {save_path}")

    return fig, ax
