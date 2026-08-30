import scanpy as sc
import numpy as np
import pandas as pd
import warnings
from tqdm import tqdm
import time

from sklearn.neighbors import KDTree

import argparse

from _common import output_directory, report_resources

parser = argparse.ArgumentParser(description="Run SpatialCCCbench Baseline 1")
parser.add_argument("--analysis_dataset", help="Analysis dataset name", required=True)
parser.add_argument("--adata_path", help="Path to adata file", required=True)
parser.add_argument("--LR_ref_path", help="Path to LR ref file", required=True)
parser.add_argument("--output_root", default="result", help="Root directory for benchmark outputs")
parser.add_argument("--n_perm", help="Number of permutations", type=int, default=1000)
parser.add_argument("--random_seed", help="Random seed for permutation", type=int, default=0)
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
LR_ligand_index_gpu = None
LR_receptor_index_gpu = None

#ensure ligand and receptor genes existing in adata.var
df = pd.read_csv(LR_database_path)
LR_ref = df[['ligand.symbol','receptor.symbol']]
LR_ref.columns = ["ligand","receptor"]

def calculate_perm_p_value(adata,index_list,spot_position,tree,original_express_df,LR_ref):
    df = original_express_df
    adata_ligand = adata[range(len(df.index)),list(LR_ref["ligand_index"])]
    adata_receptor = adata[range(len(df.index)),list(LR_ref["receptor_index"])]
    #index_list = list(range(len(adata.obs)))
    #find neighbor spots of all spots in tissue#
    df_neighbor = finding_neighbor(adata,spot_position,tree,index_list)

    #calculate pvalue based on permutation test
    with tqdm(
        total=len(index_list),
        desc=f"calculating_p_value",
        bar_format="{l_bar}{bar} [ time left: {remaining} ]"
    ) as pbar:
        result = []
        for index in index_list:
            #extract analysis spot and nearest six spots
            neighbor_cell = list(df_neighbor["index"][df_neighbor[f'Neighbor_{index}']==8])
            cell = list(df_neighbor["index"][df_neighbor[f'Neighbor_{index}']==10])
            if use_gpu:
                neighbor_index = torch.as_tensor(neighbor_cell, dtype=torch.long, device=torch_device)
                ligand_express_t=expression_gpu[cell[0], LR_ligand_index_gpu]
                receptor_express_t=expression_gpu[cell[0], LR_receptor_index_gpu]
                ligand_express_ctrl=torch.mean(expression_gpu[neighbor_index[:, None], LR_ligand_index_gpu], dim=0)
                receptor_express_ctrl=torch.mean(expression_gpu[neighbor_index[:, None], LR_receptor_index_gpu], dim=0)
            else:
                ligand_express_t=adata_ligand[cell].X.toarray().ravel()
                receptor_express_t=adata_receptor[cell].X.toarray().ravel()
                ligand_express_ctrl=np.asarray(adata_ligand[neighbor_cell].X.mean(axis=0)).ravel()
                receptor_express_ctrl=np.asarray(adata_receptor[neighbor_cell].X.mean(axis=0)).ravel()
            #get random background distribution of average of product of ligand and receptor
            background=bulid_background_distribution(ligand_express_t,receptor_express_t,ligand_express_ctrl,receptor_express_ctrl,args.n_perm)
            filter_gene = adata[index].X.toarray()>np.percentile(adata[index].X.toarray(), 25)
            ligand_list = sorted(set(LR_ref["ligand"][LR_ref["ligand"].isin(adata.var.index[list(filter_gene[0])])]))
            #get pvalue of each observed LR pairs product
            for ligand in ligand_list:
                receptor_list = sorted(set(LR_ref['receptor'][LR_ref['ligand']==ligand]))
                for receptor in receptor_list:
                    if use_gpu:
                        ligand_index = gene_to_index[ligand]
                        receptor_index = gene_to_index[receptor]
                        lig_ex_t=expression_gpu[cell[0], ligand_index]
                        lig_ex_ctrl=torch.mean(expression_gpu[neighbor_index, ligand_index])
                        rec_ex_t=expression_gpu[cell[0], receptor_index]
                        rec_ex_ctrl=torch.mean(expression_gpu[neighbor_index, receptor_index])
                        observe = torch.log2((lig_ex_t+1)/(lig_ex_ctrl+1))+torch.log2((rec_ex_ctrl+1)/(rec_ex_t+1))
                    else:
                        lig_ex_t=adata[df.index[cell],ligand].X.toarray()
                        lig_ex_ctrl=np.mean(adata[df.index[neighbor_cell],ligand].X.toarray())
                        rec_ex_t =adata[df.index[cell],receptor].X.toarray()
                        rec_ex_ctrl=np.mean(adata[df.index[neighbor_cell],receptor].X.toarray())
                        lig_ex_t = lig_ex_t.item()
                        rec_ex_t = rec_ex_t.item()
                        observe =  float(np.log2((lig_ex_t+1)/(lig_ex_ctrl+1))+np.log2((rec_ex_ctrl+1)/(rec_ex_t+1)))
                    #print(ligand,receptor)
                    if use_gpu:
                        p_value = float((torch.sum(background >= observe)/len(background)).item())
                    else:
                        p_value = np.sum(background >= observe)/len(background)
                    if p_value < 0.05:
                        result.append({
                                    'ligand_Gene':ligand,
                                    'Receiver_Gene': receptor,
                                    'p_value': p_value,
                                    'sender':df.index[cell][0]
                                        })
            pbar.update(1)
    results = pd.DataFrame(result)
    return results

def finding_neighbor(adata,spot_position,tree,index_list):
    df = pd.DataFrame({"index":range(len(adata.obs))})
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

def bulid_background_distribution(ligand_express_t,receptor_express_t,ligand_express_ctrl,receptor_express_ctrl,n_perm=1000):
    log2 = torch.log2 if use_gpu else np.log2
    ligand_statistic = log2((ligand_express_t+1)/(ligand_express_ctrl+1))
    receptor_statistic = log2((receptor_express_ctrl+1)/(receptor_express_t+1))
    background = []
    for _ in range(n_perm):
        if use_gpu:
            receptor_perm = receptor_statistic[torch.randperm(receptor_statistic.numel(), device=torch_device)]
        else:
            receptor_perm = np.random.permutation(receptor_statistic)
        perm_statistic = ligand_statistic + receptor_perm
        if use_gpu:
            background.append(perm_statistic)
        else:
            background = background+list(perm_statistic)
    if use_gpu:
        return torch.cat(background)
    return np.array(background)

#ensure ligand and receptor genes existing in adata.var
df = pd.read_csv(LR_database_path)
adata_path = adata_path

destination_folder = output_directory(args, "baseline_1")

adata = sc.read_h5ad(adata_path)
adata.var_names_make_unique()

sc.pp.filter_genes(adata,min_cells=1)
sc.pp.normalize_total(adata, inplace=True)
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
if use_gpu:
    expression_gpu = torch.as_tensor(original_express_df.values, dtype=torch.float64, device=torch_device)

LR_ref = LR_ref[LR_ref['ligand'].isin(original_express_df.columns)]
LR_ref = LR_ref[LR_ref['receptor'].isin(original_express_df.columns)]
gene_to_index = {gene: idx for idx, gene in enumerate(original_express_df.columns)}
LR_ref["ligand_index"] = LR_ref["ligand"].map(gene_to_index)
LR_ref["receptor_index"] = LR_ref["receptor"].map(gene_to_index)
if use_gpu:
    LR_ligand_index_gpu = torch.as_tensor(list(LR_ref["ligand_index"]), dtype=torch.long, device=torch_device)
    LR_receptor_index_gpu = torch.as_tensor(list(LR_ref["receptor_index"]), dtype=torch.long, device=torch_device)

#deal with warning:PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.#
warnings.filterwarnings(action='ignore', category=pd.errors.PerformanceWarning)

#extarct neccessary data and create KDTree of spatial position#
spot_position = adata.obsm["spatial"]
tree = KDTree(spot_position)
index_list = list(range(len(adata.obs)))

start = time.time()
result=calculate_perm_p_value(adata,index_list,spot_position,tree,original_express_df,LR_ref)
end = time.time()
result.to_csv(destination_folder / "result.csv")
report_resources(end - start)
