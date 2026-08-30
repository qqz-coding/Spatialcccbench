import scanpy as sc
import numpy as np
import pandas as pd
import warnings
from tqdm import tqdm
import time

from sklearn.neighbors import KDTree
from scipy.stats import ttest_ind
from scipy.stats import t as t_distribution
import argparse

from _common import output_directory, report_resources

parser = argparse.ArgumentParser(description="Run SpatialCCCbench Baseline 2")
parser.add_argument("--analysis_dataset", help="Analysis dataset name", required=True)
parser.add_argument("--adata_path", help="Path to adata file", required=True)
parser.add_argument("--LR_ref_path", help="Path to LR ref file", required=True)
parser.add_argument("--output_root", default="result", help="Root directory for benchmark outputs")
parser.add_argument("--random_seed", help="Random seed", type=int, default=0)
args = parser.parse_args()
adata_path = args.adata_path
analysis_dataset = args.analysis_dataset
LR_database_path = args.LR_ref_path

np.random.seed(args.random_seed)

def get_torch_device():
    try:
        import torch

        if torch.cuda.is_available():
            torch.manual_seed(args.random_seed)
            torch.cuda.manual_seed_all(args.random_seed)
            return torch, torch.device("cuda"), True
        warnings.warn("Torch CUDA is unavailable; falling back to CPU.")
    except Exception as exc:
        warnings.warn(f"Torch is unavailable; falling back to CPU. ({exc})")
    return None, None, False

torch, torch_device, use_gpu = get_torch_device()
print(f"compute_device: {'gpu' if use_gpu else 'cpu'}")

expression_gpu = None
gene_to_index = {}
cell_to_index = {}

def ttest_ind_auto(a, b):
    if not use_gpu:
        return ttest_ind(a, b)

    n1 = a.numel()
    n2 = b.numel()
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan

    mean1 = torch.mean(a)
    mean2 = torch.mean(b)
    var1 = torch.var(a, unbiased=True)
    var2 = torch.var(b, unbiased=True)
    df = n1 + n2 - 2
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / df
    denominator = torch.sqrt(pooled_var * (1 / n1 + 1 / n2))
    if torch.isnan(denominator) or torch.isinf(denominator) or denominator.item() == 0:
        return np.nan, np.nan
    t_stat = float(((mean1 - mean2) / denominator).item())
    p_value = float(2 * t_distribution.sf(abs(t_stat), df))
    return t_stat, p_value

def ttest_gene_groups(gene, group_a, group_b):
    if not use_gpu:
        return ttest_ind(group_a[gene], group_b[gene])

    gene_index = gene_to_index[gene]
    group_a_index = torch.as_tensor([cell_to_index[cell] for cell in group_a.index], dtype=torch.long, device=torch_device)
    group_b_index = torch.as_tensor([cell_to_index[cell] for cell in group_b.index], dtype=torch.long, device=torch_device)
    a = expression_gpu[group_a_index, gene_index]
    b = expression_gpu[group_b_index, gene_index]
    return ttest_ind_auto(a, b)

def cal_spotA_sig_LR(ligand_gene,df,original_express_df,LR_ref,tree):
    '''
    ligand_gene: list(str)             specific ligand gene symbols
    adata: anndata                     adata read by scanpy from original file
    original_express_df: dataframe     cell_gene express matrix extract from adata.X
    LR_ref: dataframe                  processed LR reference database from CellChatDB
    '''

    df_ligand = pd.DataFrame(df[ligand_gene])
    df = df_ligand.copy()
    df['index'] = range(len(df))
    index_list = list(df["index"][df.iloc[:,0]!=0])
    results = []
    with tqdm(
        total=len(index_list),
        desc=f"calculate p value of each spot and neighbor spot",
        bar_format="{l_bar}{bar} [ time left: {remaining} ]"
    ) as pbar:
        if len(index_list) > 0:
        #find neighbor spot#
            df = finding_neighbor(df,spot_position,tree,index_list)

            #extract farthest_cell as control group, might be optimized in future, cause some gene dont have farthest_cell#
            df.drop(df.columns[:2], axis=1, inplace=True)
            df_neighbor_info=df
            df['max'] = df.max(axis=1)
            df_analysis = pd.DataFrame(df['max'])
            farthest_cell = df_analysis[df_analysis['max']==0]
            if len(farthest_cell) > 0:
            #calculate p value of each spot contain specific ligand gene in three grouping way#
                results = calculate_p_value(ligand_gene, df_analysis, index_list, LR_ref, original_express_df, df_neighbor_info,results)
            else:
                results.append({'ligand_Gene':ligand_gene,
		                        'Receiver_Gene': "wide_distribution_gene",
		                        't_stat': 0,
		                        'p_value': 0,
		                        'logFC':0,
		                        'group': 0,
		                        'sender': 0,
		                        'receiver': 0
		                      })
                results = pd.DataFrame(results)
        else:
            results.append({'ligand_Gene':ligand_gene,
            'Receiver_Gene': "weak_express_gene",
            't_stat': 0,
            'p_value': 0,
            'logFC':0,
            'group': 0,
            'sender': 0,
            'receiver': 0
            })
            results = pd.DataFrame(results)
        pbar.update(1)
    return results


def calculate_p_value(ligand_gene, df_analysis, index_list, LR_ref, original_express_df, df_neighbor_info,results):
    for index in index_list:
        df_gene_cell = original_express_df
        df = df_neighbor_info
        #extract farthest cell as control group,and three experimental group grouping in three way by different gene#
        farthest_cell_group = df_analysis[df_analysis['max']==0]
        sub_farthest_cell_group = df[(df[f"Neighbor_{index}"]>=4) & (df[f"Neighbor_{index}"]<10)]
        sub_nearest_cell_group = df[(df[f"Neighbor_{index}"]>=6) & (df[f"Neighbor_{index}"]<10)]
        nearest_cell_group = df[df[f"Neighbor_{index}"]==8]
        cell_analyis = df[df[f"Neighbor_{index}"]==10]

        #apply gene_expression information to diferent group
        farthest_df = df_gene_cell.loc[farthest_cell_group.index]
        nearest_df = df_gene_cell.loc[nearest_cell_group.index]
        sub_nearest_df = df_gene_cell.loc[sub_nearest_cell_group.index]
        sub_farthest_df = df_gene_cell.loc[sub_farthest_cell_group.index]

        #find corresponding receptor gene in LR_ref
        receptor_list = sorted(set(LR_ref['receptor'][LR_ref['ligand']==ligand_gene]))

        #calculate p value of each receptor gene in each neighbor spot
        for gene in receptor_list:
            t_stat, p_value = ttest_gene_groups(gene, farthest_df, nearest_df)
            logFC = np.log2((np.mean(nearest_df[gene])+1)/(np.mean(farthest_df[gene])+1))
            results.append({
                'ligand_Gene':ligand_gene,
                'Receiver_Gene': gene,
                't_stat': t_stat,
                'p_value': p_value,
                'logFC':logFC,
                'group': "nearest vs farthest",
                'sender': ', '.join(cell_analyis.index.tolist()),
                'receiver':', '.join(nearest_cell_group.index.tolist())
            })
            t_stat, p_value = ttest_gene_groups(gene, farthest_df, sub_nearest_df)
            logFC = np.log2((np.mean(sub_nearest_df[gene])+1)/(np.mean(farthest_df[gene])+1))
            results.append({
                'ligand_Gene':ligand_gene,
                'Receiver_Gene': gene,
                't_stat': t_stat,
                'p_value': p_value,
                'logFC':logFC,
                'group': "sub_nearest vs farthest",
                'sender': ', '.join(cell_analyis.index.tolist()),
                'receiver':', '.join(sub_nearest_cell_group.index.tolist())
            })
            t_stat, p_value = ttest_gene_groups(gene, farthest_df, sub_farthest_df)
            logFC = np.log2((np.mean(sub_farthest_df[gene])+1)/(np.mean(farthest_df[gene])+1))
            results.append({
                'ligand_Gene':ligand_gene,
                'Receiver_Gene': gene,
                't_stat': t_stat,
                'p_value': p_value,
                'logFC':logFC,
                'group': "sub_farthest vs farthest",
                'sender':', '.join(cell_analyis.index.tolist()),
                'receiver':', '.join(sub_farthest_cell_group.index.tolist())
            })
    results = pd.DataFrame(results)
    return results

def finding_neighbor(df,spot_position,tree,index_list):
    for index in index_list:

        distance_1, indice_1 = tree.query([spot_position[index]], k=7)
        distance_2, indice_2 = tree.query([spot_position[index]], k=19)
        distance_3, indice_3 = tree.query([spot_position[index]], k=37)

        index_1 = indice_1[0]
        index_2 = indice_2[0]
        index_3 = indice_3[0]

        df[f'Neighbor_{index}'] = 0
        df.loc[df['index'].isin(index_3), f'Neighbor_{index}'] = 4
        df.loc[df['index'].isin(index_2), f'Neighbor_{index}'] = 6
        df.loc[df['index'].isin(index_1), f'Neighbor_{index}'] = 8
        df.loc[df['index']==index, f'Neighbor_{index}'] =10
    return df

adata = sc.read_h5ad(adata_path)
adata.var_names_make_unique()

sc.pp.filter_genes(adata,min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

adata = adata[:,~adata.var["mt"].values]
adata = adata[:,~adata.var["ribo"].values]
adata = adata[:,~adata.var["hb"].values]

cell = adata.obs.index
gene = adata.var.index
original_express_df=pd.DataFrame(adata.X.toarray(), index=cell, columns=gene)
gene_to_index = {gene_name: idx for idx, gene_name in enumerate(original_express_df.columns)}
cell_to_index = {cell_name: idx for idx, cell_name in enumerate(original_express_df.index)}
if use_gpu:
    expression_gpu = torch.as_tensor(original_express_df.values, dtype=torch.float64, device=torch_device)

df = pd.read_csv(LR_database_path)
LR_ref = df[['ligand.symbol','receptor.symbol']]
LR_ref.columns = ["ligand","receptor"]

LR_ref = LR_ref[LR_ref["ligand"].isin(adata.var.index)]
LR_ref = LR_ref[LR_ref["receptor"].isin(adata.var.index)]

#deal with warning:PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.#
warnings.filterwarnings(action='ignore', category=pd.errors.PerformanceWarning)

#extarct neccessary data#
df = original_express_df.copy()
df['index'] = range(len(df))
spot_position = adata.obsm["spatial"]
tree = KDTree(spot_position)

start = time.time()
gene_count = 0
ligand_list = set(LR_ref['ligand'])
total_count = len(set(ligand_list))
results=pd.DataFrame()

print(len(ligand_list))

for ligand_gene in sorted(ligand_list):
   result = []
   result=cal_spotA_sig_LR(ligand_gene,df,original_express_df,LR_ref,tree)
   results = pd.concat([results, result], axis=0,ignore_index=True)
   print(f"{gene_count}/{total_count}_{ligand_gene}_finished")
   gene_count = gene_count+1
end = time.time()

destination_folder = output_directory(args, "baseline_2")

results.to_csv(destination_folder / "result.csv")
report_resources(end - start)
