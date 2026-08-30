import pandas as pd
import numpy as np
import anndata as ann
import scanpy as sc

import spatialdm as sdm
from spatialdm.datasets import dataset
import spatialdm.plottings as pl

import time

from scipy.sparse import csc_matrix

import matplotlib.pyplot as plt
from tqdm import tqdm
from itertools import zip_longest
import argparse

from _common import get_cell_type_matrix, output_directory, report_resources

parser = argparse.ArgumentParser(description="Run the manuscript SpatialDM workflow")
parser.add_argument("--analysis_dataset", help="Analysis dataset name", required=True)
parser.add_argument("--adata_path", help="Path to adata file", required=True)
parser.add_argument("--LR_ref_path", help="Path to LR ref file", required=True)
parser.add_argument("--output_root", default="result", help="Root directory for benchmark outputs")
parser.add_argument("--cluster_key", default="cell_type")
parser.add_argument("--n_perms", type=int, default=1000)
parser.add_argument("--n_processes", type=int, default=1)
args = parser.parse_args()
adata_path = args.adata_path
analysis_dataset = args.analysis_dataset
LR_database_path = args.LR_ref_path

adata = sc.read_h5ad(adata_path)
#adata = sc.read_visium("./example/dataset/V1_Human_Lymph_Node")
adata.var_names_make_unique()
adata.raw = adata

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

cell_type_matrix = get_cell_type_matrix(adata, args.cluster_key)
adata.obs = cell_type_matrix.copy()
# Build the SpatialDM weight matrix.
sdm.weight_matrix(adata, l=1.2, cutoff=0.2, single_cell=False) # weight_matrix by rbf kernel

interaction = pd.read_csv(LR_database_path)
df_1= interaction[['interaction_name','ligand.symbol','receptor.symbol',"annotation"]]
df_1['receptor.symbol']=df_1['receptor.symbol'].str.split(', ')
df_2= df_1.explode('receptor.symbol').reset_index(drop=True)
df_2["lrs"]=df_2.apply(lambda row:f"{row['ligand.symbol']}_{row['receptor.symbol']}",axis=1)

df_L = pd.DataFrame({"Ligand0": list(df_2["ligand.symbol"])},index = df_2["lrs"])
df_R = pd.DataFrame({"Receptor0": list(df_2["receptor.symbol"])},index = df_2["lrs"])

seq_result = list(adata.var.index)
ligand_gene = df_L["Ligand0"].values
receptor_gene = df_R["Receptor0"].values

ligand_1 = set(seq_result).intersection(set(ligand_gene))
df=df_L[df_L["Ligand0"].isin(ligand_1)]
ligand_pair=df.index

receptor_1 = set(seq_result).intersection(set(receptor_gene))
df=df_R[df_R["Receptor0"].isin(receptor_1)]
receptor_pair = df.index

total_pair = set(receptor_pair).intersection(set(ligand_pair))
df_L = df_L[df_L.index.isin(total_pair)]
df_R = df_R[df_R.index.isin(total_pair)]
df_inter_1 = pd.DataFrame({"interaction_name": list(df_2["lrs"]),"annotation":list(df_2["annotation"])},index=df_2["lrs"])
df_inter_1 = df_inter_1[df_inter_1.index.isin(total_pair)]
geneInter = df_inter_1

geneInter["ligand"]=list(df_L["Ligand0"])
geneInter["receptor"]=list(df_R["Receptor0"])
min_cell = 3
mean = "algebra"

geneInter = geneInter.sort_values('annotation')
ligand = geneInter.ligand.values
receptor = geneInter.receptor.values
geneInter.pop('ligand')
geneInter.pop('receptor')

t = []
for i in tqdm(range(len(ligand))):
    if (len(ligand[i]) > 0) * (len(receptor[i]) > 0):
        if mean=='geometric':
            meanL = gmean(adata[:, ligand[i]].X, axis=1)
            meanR = gmean(adata[:, receptor[i]].X, axis=1)
        else:
            meanL = adata[:, ligand[i]].X.mean(axis=1)
            meanR = adata[:, receptor[i]].X.mean(axis=1)
        if (sum(meanL > 0) >= min_cell) * \
                (sum(meanR > 0) >= min_cell):
            t.append(True)
        else:
            t.append(False)
    else:
        t.append(False)

ind = geneInter[t].index
adata.uns['ligand'] = pd.DataFrame.from_records(zip_longest(*pd.Series(ligand[t]).values)).transpose()
adata.uns['ligand'].columns = ['Ligand' + str(i) for i in range(adata.uns['ligand'].shape[1])]
adata.uns['ligand'].index = ind
adata.uns['receptor'] = pd.DataFrame.from_records(zip_longest(*pd.Series(receptor[t]).values)).transpose()
adata.uns['receptor'].columns = ['Receptor' + str(i) for i in range(adata.uns['receptor'].shape[1])]
adata.uns['receptor'].index = ind
adata.uns['num_pairs'] = len(ind)
adata.uns['geneInter'] = geneInter.loc[ind]

adata.uns["mean"] = mean
df_L_1=df_L[t]
df_R_1=df_R[t]
df_inter_2=df_inter_1[t]
df_inter_2["ligand"]=df_L_1["Ligand0"]
df_inter_2["receptor"] = df_R_1["Receptor0"]
df = df_L_1
df_L_1= df_L_1[~df_L_1.index.duplicated(keep="first")]
df_R_1= df_R_1[~df_R_1.index.duplicated(keep="first")]
geneInter= geneInter[~geneInter.index.duplicated(keep="first")]
adata.uns["ligand"]=df_L_1
adata.uns["receptor"]=df_R_1
adata.uns["num_pairs"]=len(df_L_1)
adata.uns["geneInter"]=geneInter

geneInter=geneInter[geneInter.index.isin(df_L_1.index)]
adata.uns['geneInter']=geneInter

# Select globally significant ligand-receptor pairs.
start = time.time()
sdm.spatialdm_global(adata, args.n_perms, specified_ind=None, method='both', nproc=args.n_processes)     # global Moran selection
sdm.sig_pairs(adata, method='permutation', fdr=False, threshold=0.05)     # select significant pairs
print("%.3f seconds" %(time.time()-start))

# Select locally significant spots.
sdm.spatialdm_local(adata, n_perm=args.n_perms, method='both', specified_ind=None, nproc=args.n_processes)     # local spot selection
sdm.sig_spots(adata, method='permutation', fdr=False, threshold=0.05)     # significant local spots


adata.obsm['celltypes'] = adata.obs[adata.obs.columns]

def ligand_ct(adata, pair):
    ct_L = (
        adata.uns['local_stat']['local_I'][:,adata.uns['selected_spots'].index==pair] *
        adata.obsm['celltypes']
    )
    return ct_L

def receptor_ct(adata, pair):
    ct_R = (
        adata.uns['local_stat']['local_I_R'][:,adata.uns['selected_spots'].index==pair] *
        adata.obsm['celltypes']
    )
    return ct_R

min_quantile=0.5
df_dict={}


pairs=list(adata.uns["local_perm_p"].index)
if type(min_quantile) is float:
    min_quantile = np.repeat(min_quantile, len(pairs))
for i, pair in enumerate(pairs):
    Links = pd.DataFrame()
    type_interaction = adata.uns['geneInter'].loc[pair, 'annotation']
    if type_interaction == 'Secreted Signaling':
        w = adata.obsp['weight']
    else:
        w = adata.obsp['nearest_neighbors']
    print(pair)
    ct_L = ligand_ct(adata, pair)
    ct_R = receptor_ct(adata, pair)

    sparse_ct_sum = [[(csc_matrix(w).multiply(ct_L[n1].values).T.multiply(ct_R[n2].values)).sum() \
                      for n1 in ct_L.columns] for n2 in ct_R.columns]
    sparse_ct_sum = np.array(sparse_ct_sum)

    Links = pd.DataFrame({'source': np.tile(ct_L.columns, ct_R.shape[1]),
                          'target': np.repeat(ct_R.columns, ct_L.shape[1]),
                          'value': sparse_ct_sum.reshape(1, -1)[0]})
    df = Links[:]
    df_dict[f'{pair}'] = Links

df_sum=pd.DataFrame()
dict_cell_pairs = {}
for pair in df_dict.keys():
    tmp = df_dict[pair].copy()
    #tmp = tmp.loc[tmp["value"] >= 1000]
    tmp['pairs'] = pair
    df_sum = pd.concat([df_sum, tmp], axis=0, ignore_index=True)

for source in df_sum["source"].unique():
    for target in df_sum["target"].unique():
        key = f"{source}-{target}"
        if key not in dict_cell_pairs:
            dict_cell_pairs[key] = set()

for source in df_sum["source"].unique():
    for target in df_sum["target"].unique():
        df_test = df_sum[(df_sum["source"] == source) & (df_sum["target"] == target)]
        pairs = df_test["pairs"].unique()
        key = f"{source}-{target}"
        dict_cell_pairs[key].update(pairs)
dict_cell_pairs_2 = dict_cell_pairs

destination_folder = output_directory(args, "spatialDM")
df_sum.to_csv(destination_folder / "result.csv")
report_resources(time.time() - start)
