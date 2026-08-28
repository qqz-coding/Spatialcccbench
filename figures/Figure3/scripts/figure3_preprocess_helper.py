# Submission-aligned copy generated 2026-08-18
# Submission figure: Figure3
# Role: Preprocess helper required by boundary adaptation scripts

import pandas as pd
import anndata as ad
from tqdm import tqdm
import numpy as np
from scipy.sparse import csc_matrix
import operator

import os


class _AnnDataReader:
    read_h5ad = staticmethod(ad.read_h5ad)


sc = _AnnDataReader()

def extract_result(tool_list, analysis_dataset,spot_info):
    result_df_dict = {}
    for tool in tool_list:
        if tool == "Squidpy":
            df = extract_result_Squidpy(analysis_dataset)
            result_df_dict[tool] = df
            print(f"{tool} for {analysis_dataset} finished")
        if tool == "Squidpy+omnipath":
            df = extract_result_Squidpy_1(analysis_dataset)
            result_df_dict[tool] = df
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "CellAgentChat":
            df = extract_result_CellAgentChat(analysis_dataset)
            result_df_dict[tool] = df
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "SpaTalk":
            df = extract_result_SpaTalk(analysis_dataset)
            result_df_dict[tool] = df
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "SpatialDM":
            df = extract_result_SpatialDM(analysis_dataset)
            result_df_dict[tool] = df
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "stLearn":
            df = extract_result_stLearn(analysis_dataset)
            result_df_dict[tool] = df
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "stLearn_without_spotmixture":
            result_df_dict[tool] = extract_result_stLearn_without_spotmixture(analysis_dataset)
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "COMMOT":
            result_df_dict[tool] = extract_result_COMMOT(analysis_dataset)
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "Giotto":
            result_df_dict[tool] = extract_result_giotto(analysis_dataset)
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "Baseline_1":
            result_df_dict[tool] = extract_result_baseline_1(analysis_dataset,spot_info)
            print(f"{tool} for {analysis_dataset} finished")
        elif tool == "Baseline_2":
            result_df_dict[tool] = extract_result_baseline_2(analysis_dataset,spot_info)
            print(f"{tool} for {analysis_dataset} finished")
        elif os.path.exists(f"./{analysis_dataset}/custom_ccc_reuslt.csv"):
                result_df_dict[tool]=extract_custom_result(analysis_dataset,thereshold=0.05,op_func=operator.lt)
                print(f"{tool} for {analysis_dataset} finished")
    return result_df_dict

def extract_custom_result(analysis_dataset,thereshold=None,op_func=operator.lt):
    custom = pd.read_csv(f"./{analysis_dataset}/custom_ccc_reuslt.csv",index_col=0)
    mask = op_func(custom, thereshold)
    result = mask.stack()[mask.stack()].index.tolist()
    result_df = pd.DataFrame(result, columns=['LR_pairs', 'cell_pairs'])
    custom_result = result_df
    result = []
    for pairs in tqdm(list(set(custom_result["LR_pairs"]))):
        df = custom_result[custom_result["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result

def extract_result_SpatialDM(analysis_dataset):
    spatialdm = pd.read_csv(f"./result/{analysis_dataset}/spatialDM/result.csv")
    result = []
    spatialdm["cell_pairs"] = spatialdm["source"] + "-" +spatialdm["target"]
    for pairs in tqdm(list(set(spatialdm["pairs"]))):
        df = spatialdm[spatialdm["pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result

def extract_result_COMMOT(analysis_dataset):
    commot_adata= sc.read_h5ad(f"./result/{analysis_dataset}/COMMOT/result.h5ad")
    if analysis_dataset in ["merfish/adata_minus9_ori"]:
        df = pd.read_csv("./dataset/cellchatDB_mouse_symbol.csv")
    else:
        df = pd.read_csv("./dataset/interaction_CC.csv")
    df["LR_pairs"] = df['ligand.symbol'] + "-" + df['receptor.symbol']
    LR_list =list(df["LR_pairs"])
    adata = commot_adata
    extra_LR = []
    result = []
    for pairs in tqdm(LR_list):
        key = f'commot_cluster_spatial_permutation-cell_type-lymnode-{pairs}'
        if key not in adata.uns:
            extra_LR.append(key)
        else:
            tmp = adata.uns[key]["communication_pvalue"]
            filtered_values = []
            filtered_values = tmp[tmp < 0.05]
            stacked_values = filtered_values.stack()
            multi_index = stacked_values.index
            cell_pairs = ["-".join(map(str, idx)) for idx in multi_index]
            tmp={
                "LR_pairs":pairs.replace("-","_"),
                "cell_pairs":cell_pairs,
                "interaction_cell_num":len(list(cell_pairs))
            }
            result.append(tmp)
    result= pd.DataFrame(result)
    return result

def extract_result_Squidpy(analysis_dataset):
    Squidpy_result = pd.read_csv(f"./result/{analysis_dataset}/cellphoneDB/result.csv", header=None)
    df = Squidpy_result
    LR= pd.DataFrame()
    LR["source"] = df.iloc[3:,0]
    LR["target"] = df.iloc[3:,1]
    cell = pd.DataFrame()
    cell["cell_1"]= list(df.iloc[0,2:].T)
    cell["cell_2"]= list(df.iloc[1,2:].T)
    cell["cell_pairs"] = cell["cell_1"]+"-"+cell["cell_2"]

    value = df.iloc[3:,2:]
    value.columns = cell["cell_pairs"]
    value.index.name = None
    value.index = LR["source"]+"_"+LR["target"]

    for col in value.columns:
        value[col] = pd.to_numeric(value[col], errors="coerce")
    result = value[value < 0.05].stack().index.tolist()

    df = pd.DataFrame(result, columns=["LR_pairs", "cell_pairs"])
    df["LR_pairs"]=df["LR_pairs"].str.replace("-","_")
    Squidpy_result = df
    result = []
    for pairs in tqdm(list(set(Squidpy_result["LR_pairs"]))):
        df = Squidpy_result[Squidpy_result["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result

def extract_result_Squidpy_1(analysis_dataset):
    Squidpy_result = pd.read_csv(f"./result/{analysis_dataset}/Squidpy+omnipath/result.csv", header=None)
    df = Squidpy_result
    LR= pd.DataFrame()
    LR["source"] = df.iloc[3:,0]
    LR["target"] = df.iloc[3:,1]
    cell = pd.DataFrame()
    cell["cell_1"]= list(df.iloc[0,2:].T)
    cell["cell_2"]= list(df.iloc[1,2:].T)
    cell["cell_pairs"] = cell["cell_1"]+"-"+cell["cell_2"]

    value = df.iloc[3:,2:]
    value.columns = cell["cell_pairs"]
    value.index.name = None
    value.index = LR["source"]+"_"+LR["target"]

    for col in value.columns:
        value[col] = pd.to_numeric(value[col], errors="coerce")
    result = value[value < 0.05].stack().index.tolist()

    df = pd.DataFrame(result, columns=["LR_pairs", "cell_pairs"])
    df["LR_pairs"]=df["LR_pairs"].str.replace("-","_")
    Squidpy_result = df
    result = []
    for pairs in tqdm(list(set(Squidpy_result["LR_pairs"]))):
        df = Squidpy_result[Squidpy_result["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result

def  extract_result_stLearn(analysis_dataset):
    stlearn_result = sc.read_h5ad( f"./result/{analysis_dataset}/stlearn/result.h5ad")
    #stlearn_result= sc.read_h5ad( "./result/lymph_node/stlearn/adata_lr_analysis_lymph_without_spotmixtrue.h5")
    stlearn_result = stlearn_result.uns[ 'per_lr_cci_pvals_cell_type']
    result = []
    for pairs in tqdm(list(stlearn_result.keys())):
        tmp = stlearn_result[pairs][stlearn_result[pairs]<0.05]
        non_nan = tmp.notna().stack()
        if len(non_nan[non_nan].index.tolist())>0:
            pairs_list = []
            for cell_pairs in non_nan[non_nan].index.tolist():
                string = cell_pairs[0] + "-" + cell_pairs[1]
                pairs_list.append(string)
            result.append({ "LR_pairs": pairs.replace("-","_"),
                                             "cell_pairs" : pairs_list,
                                             "interaction_cell_num": len(pairs_list)
                                                               })
    result = pd.DataFrame(result)
    return result

def  extract_result_stLearn_without_spotmixture(analysis_dataset):
    stlearn_result= sc.read_h5ad( f"./result/{analysis_dataset}/stlearn/result_without_mixture.h5ad")
    stlearn_result = stlearn_result.uns[ 'per_lr_cci_pvals_cell_type']
    result = []
    for pairs in tqdm(list(stlearn_result.keys())):
        tmp = stlearn_result[pairs][stlearn_result[pairs]<0.05]
        non_nan = tmp.notna().stack()
        if len(non_nan[non_nan].index.tolist())>0:
            pairs_list = []
            for cell_pairs in non_nan[non_nan].index.tolist():
                string = cell_pairs[0] + "-" + cell_pairs[1]
                pairs_list.append(string)
            result.append({ "LR_pairs": pairs.replace("-","_"),
                                             "cell_pairs" : pairs_list,
                                             "interaction_cell_num": len(pairs_list)
                                                               })
    result = pd.DataFrame(result)
    return result

def extract_result_CellAgentChat(analysis_dataset):
    cellagentchat_result = pd.read_csv(f"./result/{analysis_dataset}/cellagentchat/result.csv",index_col=0)
    mask = cellagentchat_result < 0.05
    result = mask.stack()[mask.stack()].index.tolist()
    result_df = pd.DataFrame(result, columns=['LR_pairs', 'cell_pairs'])
    result_df[["cell1","cell2"]] = result_df["cell_pairs"].str.split('_', expand=True)
    result_df["cell1"]=result_df["cell1"].str.replace("-","_")
    result_df["cell2"]=result_df["cell2"].str.replace("-","_")
    result_df["cell_pairs"] = result_df["cell1"] + "-" + result_df["cell2"]
    result_df["LR_pairs"] = result_df["LR_pairs"].str.replace("-","_")

    cellagentchat_result = result_df
    result = []
    for pairs in tqdm(list(set(cellagentchat_result["LR_pairs"]))):
        df = cellagentchat_result[cellagentchat_result["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result

def extract_result_giotto(analysis_dataset):
    giotto_result = pd.read_csv(f"./result/{analysis_dataset}/giotto/result.csv",index_col=0)
    giotto_result ["LR_pairs"] = giotto_result["LR_comb"].str.replace("-","_")
    giotto_result ["cell_pairs"]=giotto_result["lig_cell_type"] + "-" + giotto_result["rec_cell_type"]
    result =[]
    for pairs in tqdm(list(set(giotto_result["LR_pairs"]))):
        df = giotto_result[giotto_result["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result

def extract_result_SpaTalk(analysis_dataset):
    #spatalk
    spatalk_result =  pd.read_csv(f"./result/{analysis_dataset}/spatalk/result.csv",index_col=0)
    spatalk_result["LR_pairs"] = spatalk_result["ligand"] + "_" + spatalk_result["receptor"]
    spatalk_result["cell_pairs"] = spatalk_result["celltype_sender"] + "-" + spatalk_result["celltype_receiver"]
    result = []
    for pairs in tqdm(list(set(spatalk_result["LR_pairs"]))):
        df = spatalk_result[spatalk_result["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result

def extract_result_baseline_1(analysis_dataset,spot_info):
    baseline_1 = pd.read_csv(f"./result/{analysis_dataset}/baseline_1/result.csv",index_col=0)
    df = baseline_1
    ref = pd.read_csv(spot_info,index_col=0)
    ref["neighbor"] = ref["neighbor"].apply(extract_neighbors)
    baseline1 = pd.merge(df, ref, on="sender", how="left")
    df["receiver"] = baseline1["neighbor"].astype(str)

    label_to_type = ref["cell_type"].to_dict()
    df["sender"] = replace_labels_with_types(df["sender"],label_to_type)
    df["receiver"] = replace_labels_with_types(df["receiver"],label_to_type)

    df = df[df["p_value"]<0.05]
    df["cell_pairs"] = df.apply(lambda row: generate_cell_pairs(row["sender"], row["receiver"]), axis=1)
    df["LR_pairs"] = df["ligand_Gene"]+"_"+df["Receiver_Gene"]
    baseline1 = df[["LR_pairs","cell_pairs"]]
    baseline1.loc[:, "cell_pairs"] = baseline1["cell_pairs"].str.split(",")
    expanded_df = baseline1.explode("cell_pairs").reset_index(drop=True)
    expanded_df["LR_cell"] = expanded_df["LR_pairs"]+"_"+expanded_df["cell_pairs"]
    baseline_1 = expanded_df.drop_duplicates(subset="LR_cell", keep="first")

    result = []
    for pairs in tqdm(list(set(baseline_1["LR_pairs"]))):
        df = baseline_1[baseline_1["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)

    result= pd.DataFrame(result)
    return result

def extract_result_baseline_2(analysis_dataset,spot_info):
    baseline_2 = pd.read_csv(f"./result/{analysis_dataset}/baseline_2/result.csv",index_col=0)
    df = baseline_2[(baseline_2["p_value"]<0.05)&(baseline_2["logFC"]>0)]
    df = df[df["Receiver_Gene"]!="wide_distribution_gene"]
    ref = pd.read_csv(spot_info,index_col=0)
    label_to_type = ref["cell_type"].to_dict()
    df["sender"] = replace_labels_with_types(list(df["sender"]),label_to_type)
    df["receiver"] = replace_labels_with_types(list(df["receiver"]),label_to_type)
    df_nearest = df[df["group"]=="nearest vs farthest"]

    df_nearest.loc[:,"cell_pairs"] = df.apply(lambda row: generate_cell_pairs(row["sender"], row["receiver"]), axis=1)
    df_nearest.loc[:,"LR_pairs"] = df_nearest["ligand_Gene"]+"_"+df_nearest["Receiver_Gene"]
    baseline2 = df_nearest[["LR_pairs","cell_pairs"]]
    baseline2.loc[:, "cell_pairs"] = baseline2["cell_pairs"].str.split(",")
    expanded_df = baseline2.explode("cell_pairs").reset_index(drop=True)
    expanded_df["LR_cell"] = expanded_df["LR_pairs"]+"_"+expanded_df["cell_pairs"]

    baseline_2 = expanded_df.drop_duplicates(subset="LR_cell", keep="first")
    result = []
    for pairs in tqdm(list(set(baseline_2["LR_pairs"]))):
        df = baseline_2[baseline_2["LR_pairs"] ==pairs]
        tmp={
        "LR_pairs":pairs.replace("-","_"),
        "cell_pairs":list(df["cell_pairs"]),
        "interaction_cell_num":len(list(df["cell_pairs"]))
         }
        result.append(tmp)
    result= pd.DataFrame(result)
    return result


def generate_cell_pairs(sender, receivers):

    receivers_list = receivers.split(",")

    pairs = [f"{sender}-{receiver}" for receiver in receivers_list]

    pairs = set(pairs)

    return ",".join(pairs)

def replace_labels_with_types(df_list, label_to_type):
    label_list = []
    for labels in df_list:
        labels_list = labels.split(",")
        types_list = [label_to_type.get(label.strip(), "Unknown") for label in labels_list]
        new_label = ",".join(types_list)
        label_list.append(new_label)
    return label_list

def extract_neighbors(neighbor_str):
    neighbor_str = neighbor_str.strip("[]'")
    neighbors = neighbor_str.replace("'", "").split(", ")
    return ", ".join(neighbors)
