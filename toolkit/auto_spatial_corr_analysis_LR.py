import torch.nn as nn
from torch.nn.parallel import DataParallel
import torch

from esda.moran import Moran
from esda.moran import Moran_Local
from esda.geary import Geary
from esda.geary_local import Geary_Local

import scanpy as sc
import pandas as pd
import numpy as np
from tqdm import tqdm

from scipy import sparse


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


def prepare_weights(adata, w_type="Diffusion"):
    """Prepare spatial weight matrix for Moran's I calculation"""

    adata.var_names_make_unique()
    _normalize_log1p_once(adata)

    if w_type == "Diffusion":
        sc.pp.neighbors(adata, n_neighbors=61, use_rep="spatial")

        matrix = adata.obsp['distances']

        coo = matrix.tocoo()

        new_data = np.zeros_like(coo.data)

        # Assign different weights based on distance ranges
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

        # Add diagonal elements
        n = matrix.shape[0]
        diag_rows = np.arange(n)
        diag_cols = np.arange(n)
        diag_data = np.full(n, 12)
        diag_coo = sparse.coo_matrix((diag_data, (diag_rows, diag_cols)), shape=matrix.shape)
        final_coo = new_coo + diag_coo

        new_csr = final_coo.tocsr()

        adata.obsp['weight_diffusion'] = new_csr

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

        adata.obsp['weight_contact'] = new_csr
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

class MoranCalculator(nn.Module):
    """Neural network module for Moran's I calculation, supporting DataParallel and permutation tests"""

    def __init__(self):
        super(MoranCalculator, self).__init__()

    def forward(self, ligand_batch, receptor_batch, weights):
        return self.calculate_morans_gearys_batch(ligand_batch, receptor_batch, weights)

    def calculate_morans_gearys_batch(self, ligand_batch, receptor_batch, weights):
        """Calculate Moran's I and Geary's C statistics for a batch of ligand-receptor pairs"""

        x_mean = torch.mean(ligand_batch, dim=1, keepdim=True)
        y_mean = torch.mean(receptor_batch, dim=1, keepdim=True)

        x_all = ligand_batch - x_mean
        y_all = receptor_batch - y_mean

        x_all_sq = torch.sum(x_all * x_all, dim=1)
        y_all_sq = torch.sum(y_all * y_all, dim=1)

        numerator = torch.einsum('bi,ij,bj->b', x_all, weights, y_all)
        denominator = torch.sqrt(x_all_sq) * torch.sqrt(y_all_sq)
        total_weights = torch.sum(weights)
        n = ligand_batch.size(1)

        morans_i = torch.where((denominator != 0) & (total_weights != 0),
                              (n / total_weights) * numerator / denominator,
                              torch.zeros_like(numerator))

        local_morans = x_all * (torch.einsum('bj,ij->bi', y_all, weights))

        # Standardize each gene and row-standardize W. The 0.5 factor gives
        # local cross-Geary a random-background reference value around 1.
        # Self-weights are retained for within-spot LR communication.
        eps = torch.finfo(ligand_batch.dtype).eps
        x_std = x_all / torch.sqrt(torch.mean(x_all ** 2, dim=1, keepdim=True)).clamp_min(eps)
        y_std = y_all / torch.sqrt(torch.mean(y_all ** 2, dim=1, keepdim=True)).clamp_min(eps)
        local_squared_diff = (x_std.unsqueeze(2) - y_std.unsqueeze(1)) ** 2
        row_sums = torch.sum(weights, dim=1, keepdim=True)
        row_standardized_weights = torch.where(
            row_sums != 0,
            weights / row_sums,
            torch.zeros_like(weights),
        )
        local_geary_c = 0.5 * torch.einsum(
            'bij,ij->bi', local_squared_diff, row_standardized_weights
        )

        # Use the same normalized scale globally and locally. With row-
        # standardized weights, this is the mean local cross-Geary value.
        geary_c = torch.mean(local_geary_c, dim=1)

        return morans_i, geary_c, local_morans, local_geary_c

def permutation_test(calculator, ligand_batch, receptor_batch, weights,
                    n_permutations=1000, device='cuda'):
    """
    Perform permutation test to calculate p-values

    Parameters:
        calculator: MoranCalculator instance
        ligand_batch: Ligand expression tensor
        receptor_batch: Receptor expression tensor
        weights: Spatial weight matrix
        n_permutations: Number of permutations
        device: Computation device

    Returns:
        p_value_moran: p-value for Moran's I
        p_value_geary: p-value for Geary's C
        perm_moran_stats: Moran's I distribution from permutation test
        perm_geary_stats: Geary's C distribution from permutation test
    """
    # Calculate original statistics
    with torch.no_grad():
        original_moran, original_geary, _, _ = calculator(
            ligand_batch, receptor_batch, weights
        )

    n_pairs = ligand_batch.size(0)
    n_cells = ligand_batch.size(1)

    perm_moran_stats = torch.zeros(n_permutations, n_pairs, device=device)
    perm_geary_stats = torch.zeros(n_permutations, n_pairs, device=device)

    for perm_idx in range(n_permutations):
        # Create random permutation indices to shuffle ligand and receptor
        perm_indices = torch.randperm(n_cells, device=device)

        # Shuffle ligand expression
        perm_ligand = ligand_batch[:, perm_indices]

        # Shuffle receptor expression (maintaining ligand-receptor pairing)
        perm_receptor = receptor_batch[:, perm_indices]

        with torch.no_grad():
            perm_moran, perm_geary, _, _ = calculator(
                perm_ligand, perm_receptor, weights
            )

        perm_moran_stats[perm_idx] = perm_moran
        perm_geary_stats[perm_idx] = perm_geary

    # Calculate p-values
    p_value_moran = torch.zeros(n_pairs, device=device)
    p_value_geary = torch.zeros(n_pairs, device=device)

    for i in range(n_pairs):
        # Moran's I p-value: position of observed value in permutation distribution
        moran_extreme = torch.sum(perm_moran_stats[:, i] >= original_moran[i])
        p_value_moran[i] = (moran_extreme + 1) / (n_permutations + 1)  # Avoid p=0

        # Geary's C p-value: position of observed value in permutation distribution
        geary_extreme = torch.sum(perm_geary_stats[:, i] <= original_geary[i])
        p_value_geary[i] = (geary_extreme + 1) / (n_permutations + 1)

    return p_value_moran, p_value_geary, perm_moran_stats, perm_geary_stats

def parallel_moran_calculation(
    LR_ref, adata, weight, batch_size=64, n_permutations=1000,
    gini_threshold=0,
):
    """
    Perform multi-GPU parallel calculation using DataParallel, including permutation tests

    Parameters:
        LR_ref: Ligand-receptor reference data
        adata: AnnData object containing expression data
        weight: Spatial weight matrix
        batch_size: Batch size
        n_permutations: Number of permutations for permutation test
        gini_threshold: Minimum Gini coefficient required for both genes.
            Set to None to disable Gini filtering.

    Returns:
        result_all_LR: List containing statistics and p-values
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    num_gpus = torch.cuda.device_count()
    if num_gpus > 1:
        print(f"Using {num_gpus} GPUs for parallel calculation")

    calculator = MoranCalculator()
    if num_gpus > 1:
        calculator = DataParallel(calculator)
    calculator = calculator.to(device)

    weight_tensor = torch.as_tensor(weight, dtype=torch.float32, device=device)

    # Get number of observations (cells) in adata
    n_cells = adata.shape[0]

    unique_pairs = LR_ref[['ligand', 'receptor']].drop_duplicates()
    result_all_LR = []

    # Filter out non-existent gene pairs
    valid_pairs = []
    for _, row in unique_pairs.iterrows():
        ligand, receptor = row['ligand'], row['receptor']

        # Check if genes exist
        if ligand in adata.var_names and receptor in adata.var_names:
            valid_pairs.append((ligand, receptor))

    for i in tqdm(range(0, len(valid_pairs), batch_size)):
        batch_pairs = valid_pairs[i:i + batch_size]

        # Ensure all vectors have consistent length
        ligand_data_list = []
        receptor_data_list = []
        ligand_gini_list = []
        receptor_gini_list = []
        ligands = []
        receptors = []

        for ligand, receptor in batch_pairs:
            try:
                ligand_idx = np.where(adata.var_names == ligand)[0]
                receptor_idx = np.where(adata.var_names == receptor)[0]
                if len(ligand_idx) == 0 or len(receptor_idx) == 0:
                    continue

                ligand_col = adata.X[:, ligand_idx[0]]
                receptor_col = adata.X[:, receptor_idx[0]]
                ligand_data = (ligand_col.toarray() if sparse.issparse(ligand_col)
                               else np.asarray(ligand_col)).ravel()
                receptor_data = (receptor_col.toarray() if sparse.issparse(receptor_col)
                                 else np.asarray(receptor_col)).ravel()
                if len(ligand_data) != n_cells or len(receptor_data) != n_cells:
                    continue

                ligand_gini = gini_coefficient(ligand_data)
                receptor_gini = gini_coefficient(receptor_data)
                if gini_threshold is not None and (
                    ligand_gini < gini_threshold or receptor_gini < gini_threshold
                ):
                    continue

                ligand_data_list.append(ligand_data)
                receptor_data_list.append(receptor_data)
                ligand_gini_list.append(ligand_gini)
                receptor_gini_list.append(receptor_gini)
                ligands.append(ligand)
                receptors.append(receptor)

            except Exception as e:
                print(f"Error extracting expression for gene {ligand}-{receptor}: {e}")
                continue

        if not ligand_data_list or not receptor_data_list:
            continue

        # Check if all vectors have consistent length
        lengths = [len(ld) for ld in ligand_data_list] + [len(rd) for rd in receptor_data_list]
        if len(set(lengths)) != 1:
            print(f"Inconsistent vector lengths in batch {i//batch_size}: {set(lengths)}")
            # Pad or truncate all vectors to same length
            target_length = n_cells
            ligand_data_list = [ld[:target_length] if len(ld) >= target_length
                              else np.pad(ld, (0, target_length - len(ld)), 'constant')
                              for ld in ligand_data_list]
            receptor_data_list = [rd[:target_length] if len(rd) >= target_length
                                else np.pad(rd, (0, target_length - len(rd)), 'constant')
                                for rd in receptor_data_list]

        # Prepare batch tensors
        try:
            ligand_batch = torch.stack([torch.as_tensor(ld, dtype=torch.float32)
                                      for ld in ligand_data_list]).to(device)
            receptor_batch = torch.stack([torch.as_tensor(rd, dtype=torch.float32)
                                        for rd in receptor_data_list]).to(device)
        except Exception as e:
            print(f"Error creating batch tensors: {e}")
            continue

        # Calculate actual statistics
        with torch.no_grad():
            try:
                moran_batch, geary_batch, local_moran_batch, local_geary_batch = calculator(
                    ligand_batch, receptor_batch, weight_tensor
                )
            except Exception as e:
                print(f"Error calculating statistics: {e}")
                continue

        # Perform permutation test to calculate p-values
        try:
            p_value_moran, p_value_geary, _, _ = permutation_test(
                calculator, ligand_batch, receptor_batch, weight_tensor,
                n_permutations=n_permutations, device=device
            )
        except Exception as e:
            print(f"Error performing permutation test: {e}")
            continue

        # Collect results
        for j, (ligand, receptor) in enumerate(zip(ligands, receptors)):
            result_all_LR.append({
                "ligand": ligand,
                "receptor": receptor,
                "ligand_gini": ligand_gini_list[j],
                "receptor_gini": receptor_gini_list[j],
                "global_moran": moran_batch[j].cpu().item(),
                "global_geary": geary_batch[j].cpu().item(),
                "local_moran": local_moran_batch[j].cpu().tolist(),
                "local_geary": local_geary_batch[j].cpu().tolist(),
                "p_value_moran": p_value_moran[j].cpu().item(),
                "p_value_geary": p_value_geary[j].cpu().item(),
                "n_permutations": n_permutations
            })

    return result_all_LR

def obtain_diff_LR_df(edge_strong_ligand, group, LR_ref, result_LR):
    """Create DataFrame indicating differential ligand-receptor interactions by group"""

    new_columns_dict = {}

    for ligand in tqdm(edge_strong_ligand):
        for receptor in set(LR_ref['receptor'][LR_ref['ligand'] == ligand]):
            if receptor in set(result_LR['receptor']):
                local_morani = ast.literal_eval(result_LR['local_moran'][result_LR['ligand'] == ligand].values[0])
                local_greay = ast.literal_eval(result_LR['local_geary'][result_LR['ligand'] == ligand].values[0])

                if group == "edge":
                    bool_list = [a > 0 and b >= 1.5 for a, b in zip(local_morani, local_greay)]
                elif group == "gradient":
                    bool_list = [a > 0 and b > 0 and b <= 0.5 for a, b in zip(local_morani, local_greay)]
                elif group == "normal":
                    bool_list = [a > 0 and b > 0.5 and b < 1.5 for a, b in zip(local_morani, local_greay)]
                else:
                    bool_list = [a < 0 or b == 0 for a, b in zip(local_morani, local_greay)]

                col_name = f"{ligand}_{receptor}"
                new_columns_dict[col_name] = [1 if condition else 0 for condition in bool_list]

    if new_columns_dict:
        new_columns_df = pd.DataFrame(new_columns_dict, index=adata.obs.index)

    return new_columns_df
