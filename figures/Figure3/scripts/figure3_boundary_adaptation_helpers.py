from pathlib import Path
# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure3
# Role: Figure 3 helper functions for edge equalization, gradient change, and edge signal loss

import pandas as pd
from figure3_preprocess_helper import *

import os
from sklearn.neighbors import KDTree
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import pearsonr, probplot
import scipy.stats as stats

FIGURE_DIR = Path(os.environ.get("SPATIALCCCBENCH_OUTPUT_DIR", "figure"))


def _save_figure(fig, stem):
    os.makedirs(FIGURE_DIR, exist_ok=True)
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}.svg", dpi=300, bbox_inches="tight")


def _use_arial_font():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.unicode_minus": False,
    })


def _scatter_with_ci(ax, x, y, color="#4682B4"):
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    ax.scatter(x_values, y_values, alpha=0.5, color=color, zorder=3)

    if len(x_values) < 3 or np.allclose(x_values, x_values[0]):
        return

    order = np.argsort(x_values)
    x_sorted = x_values[order]
    y_sorted = y_values[order]
    slope, intercept = np.polyfit(x_sorted, y_sorted, 1)
    fitted = slope * x_sorted + intercept
    residual = y_sorted - fitted
    dof = len(x_sorted) - 2
    if dof <= 0:
        return

    s_err = np.sqrt(np.sum(residual ** 2) / dof)
    x_grid = np.linspace(x_sorted.min(), x_sorted.max(), 200)
    y_grid = slope * x_grid + intercept
    x_mean = x_sorted.mean()
    sxx = np.sum((x_sorted - x_mean) ** 2)
    if sxx == 0:
        return

    t_value = stats.t.ppf(0.975, dof)
    ci = t_value * s_err * np.sqrt(1 / len(x_sorted) + (x_grid - x_mean) ** 2 / sxx)
    ax.fill_between(x_grid, y_grid - ci, y_grid + ci, color=color, alpha=0.15, linewidth=0, zorder=1)
    ax.plot(x_grid, y_grid, color=color, linewidth=2, alpha=0.9, zorder=2)


def analysis_edge_adoption_dropout(tool,all_result):
    _use_arial_font()
    result_noise_dict={}
    for noise in ["02","04","06","08"]:
        DLPFC = all_result["DLPFC"][tool]
        DLPFC_down = all_result[f"DLPFC_down_{noise}"][tool]
        DLPFC["cell_pairs"]=DLPFC["cell_pairs"].astype(str)
        DLPFC_down["cell_pairs"]=DLPFC_down["cell_pairs"].astype(str)
        #DLPFC=DLPFC[DLPFC['cell_pairs'].str.contains('Layer4-Layer5', case=False, na=False)]
        #DLPFC_down=DLPFC_down[DLPFC_down['cell_pairs'].str.contains('Layer4-Layer5', case=False, na=False)]
        DLPFC_down["LR_pairs"] = DLPFC_down["LR_pairs"].replace("HLA_","HLA-",regex=True)
        DLPFC["LR_pairs"] = DLPFC["LR_pairs"].replace("HLA_","HLA-",regex=True)
        DLPFC_down[['L', 'R']]=0
        DLPFC[['L', 'R']]=0
        DLPFC_down[['L', 'R']] = DLPFC_down['LR_pairs'].str.split('_', n=1, expand=True)
        DLPFC[['L', 'R']] = DLPFC['LR_pairs'].str.split('_', n=1, expand=True)
        ligand = DLPFC['L'].value_counts()
        ligand_down= DLPFC_down['L'].value_counts()
        ligand_down  =ligand_down.reindex(ligand.index)
        ligand=ligand.fillna(0)
        ligand_down=ligand_down.fillna(0)
        df = pd.DataFrame({"ligand": ligand, "ligand_down": ligand_down})
        result_noise_dict[noise]=df

    import matplotlib.pyplot as plt
    import numpy as np

    # 假设你已经有了result_noise_dict数据
    # 提取前20个基因作为x轴
    genes = list(result_noise_dict["02"]["ligand"].keys())[:20]  # x轴坐标
    values_ori = [result_noise_dict["02"]["ligand"][gene] for gene in genes]
    values_02 = [result_noise_dict["02"]["ligand_down"].get(gene, 0) for gene in genes]
    values_04 = [result_noise_dict["04"]["ligand_down"].get(gene, 0) for gene in genes]
    values_06 = [result_noise_dict["06"]["ligand_down"].get(gene, 0) for gene in genes]
    values_08 = [result_noise_dict["08"]["ligand_down"].get(gene, 0) for gene in genes]

    # 创建折线图
    plt.figure(figsize=(10, 8))


    colors = ['#08306b', '#2171b5', '#6baed6', '#bdd7e7', '#eff3ff']  # 深蓝到浅蓝的渐变

    # 绘制各条折线 - 噪声越小颜色越深
    plt.plot(genes, values_ori, label='Original', marker='o', linewidth=2.5, color='red')  # 原始数据用红色突出
    plt.plot(genes, values_02, label='20% Noise', marker='o', linewidth=2, color=colors[0])  # 最深色 - 最小噪声
    plt.plot(genes, values_04, label='40% Noise', marker='o', linewidth=2, color=colors[1])
    plt.plot(genes, values_06, label='60% Noise', marker='o', linewidth=2, color=colors[2])
    plt.plot(genes, values_08, label='80% Noise', marker='o', linewidth=2, color=colors[3])  # 最浅色 - 最大噪声

    # 添加图表元素
    plt.title(tool,
              fontsize=25, fontweight='bold', pad=20)
    plt.xlabel('Top 20 Genes', fontsize=12)
    plt.ylabel('Receptor number', fontsize=12)
    plt.legend(loc='best', fontsize=10)

    # 调整x轴标签
    plt.xticks(rotation=45, ha='right')  # 旋转标签以避免重叠

    # 添加网格
    plt.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    _save_figure(plt.gcf(), f"{tool}_edge_dropout_plot")

    return result_noise_dict

def analysis_edge_adoption_equal(tool,all_result):
    _use_arial_font()
    DLPFC = all_result["DLPFC"][tool]
    DLPFC_equal = all_result["DLPFC_equal"][tool]
    DLPFC_equal[['L', 'R']]=0
    DLPFC[['L', 'R']]=0
    DLPFC_equal["LR_pairs"] = DLPFC_equal["LR_pairs"].replace("HLA_","HLA-",regex=True)
    DLPFC["LR_pairs"] = DLPFC["LR_pairs"].replace("HLA_","HLA-",regex=True)
    DLPFC_equal[['L', 'R']] = DLPFC_equal['LR_pairs'].str.split('_', n=1, expand=True)
    DLPFC[['L', 'R']] = DLPFC['LR_pairs'].str.split('_', n=1, expand=True)
    intersect = set(DLPFC["L"]).intersection(set(DLPFC_equal["L"]))
    DLPFC = DLPFC[DLPFC["L"].isin(intersect)]
    DLPFC_equal = DLPFC_equal[DLPFC_equal["L"].isin(intersect)]
    ligand = DLPFC['L'].value_counts()
    ligand_equal = DLPFC_equal['L'].value_counts()
    ligand_equal  =ligand_equal.reindex(ligand.index)
    ligand=ligand.fillna(0)
    ligand_equal=ligand_equal.fillna(0)
    scipy_corr, p_value = pearsonr( ligand.values,ligand_equal.values)
    plt.style.use('bmh')
    _use_arial_font()

    ax = plt.gca()
    _scatter_with_ci(ax, ligand.values, ligand_equal.values, color="#4682B4")
    plt.xlabel("Original ligand signal", fontsize=12)
    plt.ylabel("Ligand signal under edge equalization", fontsize=12)
    plt.title(f"{tool}", fontweight='bold',fontsize=25,loc="left")

    plt.text(0.05, 0.95,
             f"PCC={scipy_corr:.4e} ",
             transform=plt.gca().transAxes,
             fontsize=12,
             fontweight='bold',
             verticalalignment='top')

    plt.tight_layout()  # 自动调整布局
    _save_figure(plt.gcf(), f"{tool}_edge_equalization_plot")

def analysis_edge_adoption_gradient_change(tool,all_result):
    _use_arial_font()
    DLPFC = all_result["DLPFC"][tool]
    DLPFC_equal = all_result["DLPFC_REVERS"][tool]
    DLPFC_equal[['L', 'R']]=0
    DLPFC[['L', 'R']]=0
    DLPFC_equal["LR_pairs"] = DLPFC_equal["LR_pairs"].replace("HLA_","HLA-",regex=True)
    DLPFC["LR_pairs"] = DLPFC["LR_pairs"].replace("HLA_","HLA-",regex=True)
    DLPFC_equal[['L', 'R']] = DLPFC_equal['LR_pairs'].str.split('_', n=1, expand=True)
    DLPFC[['L', 'R']] = DLPFC['LR_pairs'].str.split('_', n=1, expand=True)
    intersect = set(DLPFC["L"]).intersection(set(DLPFC_equal["L"]))
    DLPFC = DLPFC[DLPFC["L"].isin(intersect)]
    DLPFC_equal = DLPFC_equal[DLPFC_equal["L"].isin(intersect)]
    ligand = DLPFC['L'].value_counts()
    ligand_equal = DLPFC_equal['L'].value_counts()
    ligand_equal  =ligand_equal.reindex(ligand.index)
    ligand=ligand.fillna(0)
    ligand_equal=ligand_equal.fillna(0)
    scipy_corr, p_value = pearsonr( ligand.values,ligand_equal.values)
    plt.style.use('bmh')
    _use_arial_font()

    ax = plt.gca()
    _scatter_with_ci(ax, ligand.values, ligand_equal.values, color="#D62728")
    plt.xlabel("Original ligand signal", fontsize=12)
    plt.ylabel("Ligand signal under gradient change", fontsize=12)
    plt.title(f"{tool}", fontweight='bold',fontsize=25,loc="left")

    # 在左上角添加PCC值 (轴坐标系定位)
    plt.text(0.05, 0.95,             # x,y位置 (5%左, 95%上)
             f"PCC={scipy_corr:.4e} ",# 显示文本
             transform=plt.gca().transAxes,  # 使用轴坐标系
             fontsize=12,             # 字体大小
             fontweight='bold',
             verticalalignment='top')

    plt.tight_layout()  # 自动调整布局
    _save_figure(plt.gcf(), f"{tool}_gradient_change_plot")

def deep_analysis(all_result,tool,analysis_target="DLPFC_REVERS"):
    _use_arial_font()
    DLPFC = all_result["DLPFC"][tool]
    DLPFC_equal = all_result[analysis_target][tool]

    DLPFC_equal[['L', 'R']] = 0
    DLPFC[['L', 'R']] = 0
    DLPFC_equal["LR_pairs"] = DLPFC_equal["LR_pairs"].replace("HLA_", "HLA-", regex=True)
    DLPFC["LR_pairs"] = DLPFC["LR_pairs"].replace("HLA_", "HLA-", regex=True)
    DLPFC_equal[['L', 'R']] = DLPFC_equal['LR_pairs'].str.split('_', n=1, expand=True)
    DLPFC[['L', 'R']] = DLPFC['LR_pairs'].str.split('_', n=1, expand=True)

    intersect = set(DLPFC["L"]).intersection(set(DLPFC_equal["L"]))
    DLPFC = DLPFC[DLPFC["L"].isin(intersect)]
    DLPFC_equal = DLPFC_equal[DLPFC_equal["L"].isin(intersect)]

    ligand = DLPFC['L'].value_counts()
    ligand_equal = DLPFC_equal['L'].value_counts()
    ligand_equal = ligand_equal.reindex(ligand.index)
    ligand = ligand.fillna(0)
    ligand_equal = ligand_equal.fillna(0)

    scipy_corr, p_value = pearsonr(ligand.values, ligand_equal.values)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    _scatter_with_ci(ax1, ligand.values, ligand_equal.values, color="#4682B4")
    ax1.set_xlabel("Original ligand signal", fontsize=12)
    ax1.set_ylabel("Processed ligand signal", fontsize=12)
    ax1.set_title(f"{tool} - Ligand Counts Scatter", fontweight='bold')
    ax1.text(0.05, 0.95, f"PCC={scipy_corr:.4e}",
             transform=ax1.transAxes, fontsize=12, fontweight='bold',
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    ax1.grid(True, alpha=0.3)

    data1 = ligand.values
    data2 = ligand_equal.values

    ks_statistic, ks_pvalue = stats.ks_2samp(data1, data2)


    sorted_data1 = np.sort(data1)
    sorted_data2 = np.sort(data2)

    n = min(len(sorted_data1), len(sorted_data2))
    quantiles = np.linspace(0.01, 0.99, n)

    quantiles_data1 = np.quantile(sorted_data1, quantiles)
    quantiles_data2 = np.quantile(sorted_data2, quantiles)

    ax2.scatter(quantiles_data1, quantiles_data2, alpha=0.6, color='coral')

    min_val = min(quantiles_data1.min(), quantiles_data2.min())
    max_val = max(quantiles_data1.max(), quantiles_data2.max())
    ax2.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label='y=x reference')

    if len(quantiles_data1) > 1:
        z = np.polyfit(quantiles_data1, quantiles_data2, 1)
        p = np.poly1d(z)
        ax2.plot(quantiles_data1, p(quantiles_data1), 'g-', alpha=0.8, linewidth=2,
                label=f'Fit: y={z[0]:.3f}x+{z[1]:.3f}')

    ax2.set_xlabel('Original Ligand Count Quantiles', fontsize=12)
    ax2.set_ylabel('Gradient Ligand Count Quantiles', fontsize=12)
    ax2.set_title(f'{tool} - Q-Q Plot (Distribution Comparison)\n({analysis_target})', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    ax2.text(0.05, 0.95, f'KS statistic = {ks_statistic:.4f}\nKS p-value = {ks_pvalue:.4e}',
            transform=ax2.transAxes, fontsize=10, fontweight='bold',
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.8))

    plt.tight_layout()
    _save_figure(fig, f"{tool}_comparison_{analysis_target}")

    fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(15, 6))

    stats.probplot(data1, dist="norm", plot=ax3)
    ax3.set_title(f'{tool} - Original Data Normality Q-Q Plot\n({analysis_target})', fontweight='bold')

    stats.probplot(data2, dist="norm", plot=ax4)
    ax4.set_title(f'{tool} - Data Normality Q-Q Plot\n({analysis_target})', fontweight='bold')

    plt.tight_layout()
    _save_figure(plt.gcf(), f"{tool}_normality_qq_{analysis_target}")

    if ks_pvalue > 0.05:
        print("Distribution dont exist difference (p > 0.05)")
    else:
        print("Distribution exist difference (p <= 0.05)")
    return ks_statistic, ks_pvalue
