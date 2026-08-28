import pandas as pd
import itertools
from functools import reduce

def obtain_LR_cell_list(analysis_dict):

    LR_cell_list = {}
    for tool in list(analysis_dict.keys()):
        df =analysis_dict[tool]
        result_list = []
        for index, row in df.iterrows():
            lr_pair = row["LR_pairs"]
            cell_pairs_list = row["cell_pairs"]

            # 组合 LR_pairs 和每个 cell_pair
            for cell_pair in cell_pairs_list:
                combined_str = f"{lr_pair}--{cell_pair}"
                result_list.append(combined_str)
        LR_cell_list[tool] = result_list
        print(len(LR_cell_list[tool]))
        print(f"{tool}_finished")

    return LR_cell_list


def  obtain_all_result(LR_cell_list):
    results = []
    tools = list(LR_cell_list.keys())

    for n in range(1, len(tools) + 1):
        for combination in itertools.combinations(tools, n):
            intersection_result = set(LR_cell_list[combination[0]])
            for tool in combination[1:]:
                intersection_result = intersection_result.intersection(LR_cell_list[tool])
            results.append({
                "Combination": ", ".join(combination),
                "Intersection": intersection_result,
                "Count": len(intersection_result),
                "tool_count": n
            })
    result_df = pd.DataFrame(results)

    return result_df

def  get_all_LR_cellpair_precision(df,ground_truth_level=3):
    intersections = set()

    for item in df.loc[df["tool_count"]==ground_truth_level, 'Intersection']:
        intersections.update(item)

    unique_intersections = list(intersections)

    all_sets = df['Intersection'].tolist()

    merged_set = reduce(lambda x, y: x | y, all_sets)

    results = []

    for tool in list(df[df["tool_count"] == 1]["Combination"]):
        tool_result = df.loc[df["Combination"] == tool, "Intersection"].values[0]
        if len(tool_result) ==0:
            overlap = len(tool_result)/len(merged_set)
            precision = 0
            recall = len(tool_result.intersection(unique_intersections)) / len(unique_intersections)
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
            results.append({
                "Tool": tool,
                "overlap":overlap,
                "Precision": 0,
                "Recall": recall,
                "F1 Score": f1
            })

        else:
            overlap = len(tool_result)/len(merged_set)
            precision = len(tool_result.intersection(unique_intersections)) / len(tool_result)
            recall = len(tool_result.intersection(unique_intersections)) / len(unique_intersections)
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

            results.append({
                "Tool": tool,
                "overlap":overlap,
                "Precision": precision,
                "Recall": recall,
                "F1 Score": f1
            })
    metrics_df = pd.DataFrame(results)

    return metrics_df

def get_LR_from_per_cellpair_precision(df,ground_truth_level=3):

    results = []

    intersections = set()
    for item in df.loc[df["tool_count"]>=ground_truth_level, 'Intersection']:
        intersections.update(item)

    unique_intersections = list(intersections)
    all_sets = df['Intersection'].tolist()
    merged_set = reduce(lambda x, y: x | y, all_sets)
    split_data = [item.split('--') for item in unique_intersections]
    t_LR_CELL= pd.DataFrame(split_data, columns=['Ligand_Receptor', 'Cell_Interaction'])

    for tool in set(df[df["tool_count"] == 1]["Combination"]):
        tool_result = df.loc[df["Combination"] == tool, "Intersection"].values[0]
        split_data = [item.split('--') for item in tool_result]
        tool_LR_CELL= pd.DataFrame(split_data, columns=['Ligand_Receptor', 'Cell_Interaction'])
        #the LR precision in high confidence cell relation
        for cell_pairs in set(t_LR_CELL["Cell_Interaction"]):
            obs = set(tool_LR_CELL["Ligand_Receptor"][tool_LR_CELL['Cell_Interaction']==cell_pairs])
            grt = set(t_LR_CELL["Ligand_Receptor"][t_LR_CELL['Cell_Interaction']==cell_pairs])

            precision = len(obs.intersection(grt)) / len(obs) if len(obs)>0 else 0.0
            recall = len(obs.intersection(grt)) / len(grt)
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
            results.append({
                "Tool": tool,
                "cell_pairs":cell_pairs,
                "Precision": precision,
                "Recall": recall,
                "F1 Score": f1
            })
    metrics_df = pd.DataFrame(results)
    return metrics_df
