import os
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

from sklearn.neighbors import KDTree

from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix

import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from scipy import sparse

from esda.moran import Moran
from esda.moran import Moran_Local
from esda.geary import Geary
from esda.geary_local import Geary_Local

import scanpy as sc
import pandas as pd
import numpy as np
from tqdm import tqdm
from libpysal.weights import full2W


def _normalize_log1p_once(adata):
    """Normalize and log-transform ``adata.X`` at most once."""
    # Scanpy records successful log1p preprocessing in ``adata.uns``. This
    # also protects AnnData objects that were log-transformed before loading.
    if "log1p" in adata.uns:
        return

    if adata.raw is None:
        adata.raw = adata
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)


def prepare_weights(adata,w_type="Diffusion"):
    adata.var_names_make_unique()
    _normalize_log1p_once(adata)
    if w_type=="Diffusion":
        sc.pp.neighbors(adata, n_neighbors=61, use_rep="spatial")

        matrix = adata.obsp['distances']

        coo = matrix.tocoo()

        new_data = np.zeros_like(coo.data)

        mask_range1 = (coo.data >= 100) & (coo.data <= 200)
        mask_range2 = (coo.data >= 200) & (coo.data <= 300)
        mask_range3 = (coo.data >= 300) & (coo.data <= 400)
        mask_range4 = (coo.data >= 400) & (coo.data <= 500)
        mask_range5 = (coo.data >= 500) & (coo.data <= 600)

        new_data[mask_range1] = 10.0
        new_data[mask_range2] = 8.0
        new_data[mask_range3] = 6.0
        new_data[mask_range4] = 4.0
        new_data[mask_range5] = 2.0

        new_coo = sparse.coo_matrix((new_data, (coo.row, coo.col)), shape=matrix.shape)

        n = matrix.shape[0]
        diag_rows = np.arange(n)
        diag_cols = np.arange(n)
        diag_data = np.full(n, 12)
        diag_coo = sparse.coo_matrix((diag_data, (diag_rows, diag_cols)), shape=matrix.shape)
        final_coo = new_coo + diag_coo

        new_csr = final_coo.tocsr()

        adata.obsp['weight_diffusion']=new_csr

        weight = adata.obsp["weight_diffusion"].toarray()

    else:
        sc.pp.neighbors(adata, n_neighbors=7, use_rep="spatial")
        matrix = adata.obsp['distances']

        coo = matrix.tocoo()

        new_data = np.zeros_like(coo.data)

        mask_range1 = (coo.data >= 100) & (coo.data <= 200)

        new_data[mask_range1] = 1.0

        new_coo = sparse.coo_matrix((new_data, (coo.row, coo.col)), shape=matrix.shape)

        n = matrix.shape[0]
        diag_rows = np.arange(n)
        diag_cols = np.arange(n)
        diag_data = np.full(n, 1)
        diag_coo = sparse.coo_matrix((diag_data, (diag_rows, diag_cols)), shape=matrix.shape)
        final_coo = new_coo + diag_coo

        new_csr = final_coo.tocsr()

        adata.obsp['weight_contact']=new_csr
        weight = adata.obsp['weight_contact'].toarray()
    return weight

def gini_coefficient(values):
    """Calculate the Gini coefficient for finite non-negative values."""
    x = np.asarray(values, dtype=np.float64).ravel()

    if x.size == 0 or not np.all(np.isfinite(x)):
        return np.nan
    if np.any(x < 0):
        raise ValueError("Gini coefficient requires non-negative values.")

    total = x.sum()
    if total == 0:
        return 0.0

    x = np.sort(x)
    n = x.size
    ranks = np.arange(1, n + 1, dtype=np.float64)
    return np.sum((2 * ranks - n - 1) * x) / (n * total)

def cal_ligand_auto_corr(adata,weights,LR_ref,w_type):
    result=[]
    if w_type=="contact":
        analysis_gene=set(LR_ref["ligand"][LR_ref['annotation'].isin(['Cell-Cell Contact','ECM-Receptor'])])
    else:
        analysis_gene=set(LR_ref["ligand"][LR_ref['annotation'].isin(['Non-protein Signaling','Secreted Signaling'])])
    for ligand in tqdm(analysis_gene):
        X = adata[:,adata.var_names==ligand].X
        data =  X.toarray().flatten().tolist()
        data = np.array(data)
        if len(set(data.tolist()))>1:
            gini_index = gini_coefficient(data)
            gearys_C = Geary(data,weights)
            morani = Moran(data, weights)
            moran_local = Moran_Local(data, weights)
            gearys_local= Geary_Local(connectivity=weights,labels=True).fit(data)
            local_geary = 0.5 * np.asarray(gearys_local.localG)
            result.append({
            "ligand":ligand,
            "gini_index":gini_index,
            "morani":morani.I,
            "morani_p":morani.p_sim,
            "local_morani":moran_local.Is,
            "local_morani_z":moran_local.z_sim,
            "local_morani_p":moran_local.p_sim,
            "gearys_c":gearys_C.C,
            "local_gearysC":local_geary,
            "local_gearysC_p":gearys_local.p_sim,
           # "local_gearysC_z":gearys_local.z_sim
            })
    result = pd.DataFrame(result)
    return result

def cal_receptor_auto_corr(adata,weights,LR_ref,w_type):
    result=[]
    if w_type=="contact":
        analysis_gene=set(LR_ref["receptor"][LR_ref['annotation'].isin(['Cell-Cell Contact','ECM-Receptor'])])
    else:
        analysis_gene=set(LR_ref["receptor"][LR_ref['annotation'].isin(['Non-protein Signaling','Secreted Signaling'])])
    for receptor in tqdm(analysis_gene):
        X = adata[:,adata.var_names==receptor].X
        data =  X.toarray().flatten().tolist()
        data = np.array(data)
        if len(set(data.tolist()))>1:
            gini_index = gini_coefficient(data)
            gearys_C = Geary(data,weights)
            morani = Moran(data, weights)
            moran_local = Moran_Local(data, weights)
            gearys_local= Geary_Local(connectivity=weights,labels=True).fit(data)
            local_geary = 0.5 * np.asarray(gearys_local.localG)
            result.append({
            "receptor":receptor,
            "gini_index":gini_index,
            "morani":morani.I,
            "morani_p":morani.p_sim,
            "local_morani":moran_local.Is,
            "local_morani_z":moran_local.z_sim,
            "local_morani_p":moran_local.p_sim,
            "gearys_c":gearys_C.C,
            "local_gearysC":local_geary,
            "local_gearysC_p":gearys_local.p_sim,
           # "local_gearysC_z":gearys_local.z_sim
            })

    result = pd.DataFrame(result)
    return result
def obtain_diff_LR_df_ligand(adata, edge_strong_ligand, group, LR_ref, result_ligand, result_receptor):
    new_columns_dict = {}
    new_columns_df = pd.DataFrame()

    for ligand in tqdm(edge_strong_ligand):
        for receptor in set(LR_ref['receptor'][LR_ref['ligand'] == ligand]):
            if receptor in set(result_receptor['receptor']):
                global_moran_R = result_receptor['morani'][result_receptor['receptor'] == receptor].values[0]
                global_moran_R_p = result_receptor['morani_p'][result_receptor['receptor'] == receptor].values[0]
                gini_R = result_receptor['gini_index'][result_receptor['receptor'] == receptor].values[0]

                local_morani_L = list(result_ligand['local_morani'][result_ligand['ligand'] == ligand].values[0])
                local_greay_L = list(result_ligand['local_gearysC'][result_ligand['ligand'] == ligand].values[0])
                local_morani_R = list(result_receptor['local_morani'][result_receptor['receptor'] == receptor].values[0])
                local_greay_R = list(result_receptor['local_gearysC'][result_receptor['receptor'] == receptor].values[0])

                if group == "edge":
                    bool_list = [(a > 0 and b >= 1.5) or (c > 0 and d >= 1.5)
                                   for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]

                elif group == "gradient":
                    bool_list = [(a > 0 and b > 0 and b <= 0.5) or (c >= 0 and d > 0 and d <= 0.5)
                                     for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]
                elif group == "normal":
                    bool_list = [(a > 0 and b > 0.5 and b < 1.5) or (c > 0 and d > 0.5 and d < 1.5)
                                    for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]

                elif group == "weak":
                    bool_list = [a <= 0 or b == 0 or c < 0 or d == 0
                                for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]

                else:
                    continue

                col_name = f"{ligand}_{receptor}"
                new_columns_dict[col_name] = [1 if condition else 0 for condition in bool_list]

    if new_columns_dict:
        new_columns_df = pd.DataFrame(new_columns_dict, index=adata.obs.index)

    return new_columns_df

def obtain_diff_LR_df_receptor(adata, edge_strong_receptor, group, LR_ref, result_ligand, result_receptor):
    new_columns_dict = {}
    new_columns_df = pd.DataFrame()

    for receptor in tqdm(edge_strong_receptor):
        for ligand in set(LR_ref['ligand'][LR_ref['receptor'] == receptor]):
            if ligand in set(result_ligand['ligand']):
                local_morani_L = list(result_ligand['local_morani'][result_ligand['ligand'] == ligand].values[0])
                local_greay_L = list(result_ligand['local_gearysC'][result_ligand['ligand'] == ligand].values[0])
                local_morani_R = list(result_receptor['local_morani'][result_receptor['receptor'] == receptor].values[0])
                local_greay_R = list(result_receptor['local_gearysC'][result_receptor['receptor'] == receptor].values[0])

                if group == "edge":
                    bool_list = [(a > 0 and b >= 1.5) or (c > 0 and d >= 1.5)
                                   for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]

                elif group == "gradient":
                    bool_list = [(a > 0 and b > 0 and b <= 0.5) or (c >= 0 and d > 0 and d <= 0.5)
                                     for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]
                elif group == "normal":
                    bool_list = [(a > 0 and b > 0.5 and b < 1.5) or (c > 0 and d > 0.5 and d < 1.5)
                                    for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]

                elif group == "weak":
                    bool_list = [a <= 0 or b == 0 or c < 0 or d == 0
                                for a, b, c, d in zip(local_morani_L, local_greay_L, local_morani_R, local_greay_R)]

                else:
                    continue

                col_name = f"{ligand}_{receptor}"
                new_columns_dict[col_name] = [1 if condition else 0 for condition in bool_list]

    if new_columns_dict:
        new_columns_df = pd.DataFrame(new_columns_dict, index=adata.obs.index)

    return new_columns_df

def merge_dataframes_optimized(df1, df2):

    # 确保索引一致
    if not df1.index.equals(df2.index):
        common_index = df1.index.intersection(df2.index)
        df1 = df1.loc[common_index]
        df2 = df2.loc[common_index]

    # 分析列的情况
    common_columns = df1.columns.intersection(df2.columns)
    df1_unique_columns = df1.columns.difference(df2.columns)
    df2_unique_columns = df2.columns.difference(df1.columns)

    # 创建结果数据框
    result_df = df2.copy()

    # 处理共有列 - 使用向量化操作
    for col in common_columns:
        # 创建覆盖掩码
        mask = (df1[col] == 1) & (df2[col] == 0)
        # 应用覆盖
        result_df.loc[mask, col] = 1

    # 添加df1独有的列 - 使用一次性操作
    if len(df1_unique_columns) > 0:
        result_df = pd.concat([result_df, df1[df1_unique_columns]], axis=1)

    return result_df

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import PchipInterpolator, interp1d
import warnings
warnings.filterwarnings('ignore')

# 设置美观的matplotlib风格
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

# 颜色配置
# 默认颜色列表，如果用户没有自定义颜色，将使用这个列表
DEFAULT_COLORS = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F']

def plot_cdf_comparison(df, variable_names, variable_labels,
                       colors=None, line_styles=None, line_widths=None,
                       alpha_values=None, save_path='clean_cdf_curves_only.png',
                       title='Cumulative Distribution Functions',
                       xlabel='Count Value', ylabel='Cumulative Probability',
                       figsize=(12, 8), dpi=300, transparent=True):
    """
    绘制自定义线条颜色的CDF对比图

    参数:
    ----------
    df : pandas DataFrame
        包含数据的数据框
    variable_names : list
        要绘制的变量名列表
    variable_labels : list
        对应变量的标签列表（用于图例）
    colors : list, optional
        线条颜色列表，长度应与variable_names相同
        默认为None，将使用DEFAULT_COLORS
    line_styles : list, optional
        线条样式列表，如['-', '--', '-.', ':']，长度应与variable_names相同
    line_widths : list, optional
        线条宽度列表，长度应与variable_names相同
    alpha_values : list, optional
        透明度列表，长度应与variable_names相同
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
        # 如果颜色数量不足，循环使用提供的颜色
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

    # 绘制每条CDF曲线
    for idx, (col_name, label) in enumerate(zip(variable_names, variable_labels)):
        values = df[col_name]

        # 排序值
        sorted_values = np.sort(values)

        # 计算累积概率
        y = np.arange(1, len(sorted_values) + 1) / len(sorted_values)

        # 获取唯一值和对应的最大累积概率
        unique_values = []
        unique_y = []

        for i, val in enumerate(sorted_values):
            if i == 0 or val != sorted_values[i-1]:
                unique_values.append(val)
                unique_y.append(y[i])
            else:
                unique_y[-1] = y[i]

        unique_values = np.array(unique_values)
        unique_y = np.array(unique_y)

        # 绘制平滑的CDF曲线
        if len(unique_values) > 1:
            xs_smooth = np.linspace(unique_values.min(), unique_values.max(), 5000)

            if len(unique_values) >= 4:
                pchip = PchipInterpolator(unique_values, unique_y)
                ys_smooth = pchip(xs_smooth)
            else:
                f_linear = interp1d(unique_values, unique_y, kind='linear', fill_value='extrapolate')
                ys_smooth = f_linear(xs_smooth)

            # 使用自定义的颜色、线型和宽度
            ax.plot(xs_smooth, ys_smooth,
                   linewidth=line_widths[idx],
                   linestyle=line_styles[idx],
                   label=label,
                   color=colors[idx],
                   alpha=alpha_values[idx])

    # 设置图表属性
    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=False, loc='lower right')

    # 美化坐标轴
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

    # 设置y轴范围
    ax.set_ylim(0, 1.05)

    plt.tight_layout()

    # 保存图片
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=transparent)
    print(f"CDF曲线图已保存为: {save_path}")

    return fig, ax
