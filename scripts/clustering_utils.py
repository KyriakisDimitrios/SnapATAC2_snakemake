import os
import sys
import logging
import pandas as pd
import snapatac2 as snap
from dr_utils import load_data_to_memory

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)


def _atac_clustering_block(adata, n_features, leiden_res, batch_key, level_prefix, algorithm):
    """Runs ATAC dimensionality reduction with dynamic batch correction."""
    random_state = 20120224
    logging.info(f"[{level_prefix}] Selecting top {n_features} features...")
    snap.pp.select_features(adata, n_features=n_features)
    logging.info(f"[{level_prefix}] Running Spectral Embedding...")
    snap.tl.spectral(adata,
                     features='selected',  # Uses the Rule 3 mask
                     random_state=random_state,
                     distance_metric='cosine',
                     weighted_by_sd=True,
                     inplace=True)

    # --- Dynamic Batch Correction Switch ---
    if algorithm.lower() == 'harmony':
        logging.info(f"[{level_prefix}] Running Harmony Batch Correction on '{batch_key}'...")
        snap.pp.harmony(adata, batch=batch_key)
        rep_use = 'X_spectral_harmony'  # <-- Fixed Key
    elif algorithm.lower() == 'mnc':
        logging.info(f"[{level_prefix}] Running MNC Batch Correction on '{batch_key}'...")
        snap.pp.mnc_correct(adata, batch=batch_key)  # <-- Fixed API call
        rep_use = 'X_spectral_mnn'  # <-- Fixed Key
    else:
        logging.error(f"Unknown batch algorithm '{algorithm}'. Must be 'harmony' or 'mnc'.")
        sys.exit(1)

    logging.info(f"[{level_prefix}] Building kNN graph and UMAP using {rep_use}...")
    snap.pp.knn(adata, use_rep=rep_use)
    snap.tl.umap(adata, use_rep=rep_use)

    logging.info(f"[{level_prefix}] Running Leiden Clustering (res={leiden_res})...")
    snap.tl.leiden(adata, resolution=leiden_res)

    # Store labels with the specific level prefix
    adata.obs[f'{level_prefix}_cluster'] = adata.obs['leiden']
    del adata.obs['leiden']

    # Store UMAP coordinates directly in obs
    adata.obs[f'{level_prefix}_umap_1'] = adata.obsm['X_umap'][:, 0]
    adata.obs[f'{level_prefix}_umap_2'] = adata.obsm['X_umap'][:, 1]

    return adata


def iterative_atac_clustering(input_h5ad, output_csv, batch_key, algorithm,
                              l1_features=100000, l1_res=0.8,
                              l2_features=50000, l2_res=0.2):
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # --- LEVEL 1: Global Clustering ---
    adata = load_data_to_memory(input_h5ad)

    logging.info(f"Reading virtual AnnDataSet: {input_h5ad}")

    logging.info(f"=== Starting Level 1 (Global) Clustering using {algorithm.upper()} ===")
    _atac_clustering_block(
        adata=adata, n_features=l1_features, leiden_res=l1_res,
        batch_key=batch_key, level_prefix='L1', algorithm=algorithm
    )

    # --- LEVEL 2: Iterative Sub-Clustering ---
    logging.info("=== Starting Level 2 (Sub-cluster) Iteration ===")

    l1_clusters = adata.obs['L1_cluster'].unique()
    all_subcluster_annotations = []

    for cluster_id in l1_clusters:
        logging.info(f"Processing sub-cluster: L1_{cluster_id}")

        subset_mask = adata.obs['L1_cluster'] == cluster_id
        sub_adata = adata[subset_mask].copy()

        _atac_clustering_block(
            adata=sub_adata, n_features=l2_features, leiden_res=l2_res,
            batch_key=batch_key, level_prefix='L2', algorithm=algorithm
        )

        # Prefix the sub-clusters so they don't overlap
        sub_adata.obs['L2_cluster'] = f"{cluster_id}_" + sub_adata.obs['L2_cluster'].astype(str)

        keep_cols = ['L1_cluster', 'L1_umap_1', 'L1_umap_2',
                     'L2_cluster', 'L2_umap_1', 'L2_umap_2']

        final_cols = [c for c in sub_adata.obs.columns if c not in keep_cols] + keep_cols
        all_subcluster_annotations.append(sub_adata.obs[final_cols].copy())

        del sub_adata
        logging.info(f"Completed sub-cluster: L1_{cluster_id}")

    # --- Compile and Save ---
    logging.info("=== Compiling Final Annotations ===")
    master_obs_df = pd.concat(all_subcluster_annotations)
    master_obs_df.to_csv(output_csv, index_label='barcode')

    logging.info(f"Successfully saved {algorithm.upper()} iterative annotations to {output_csv}")

def run_l2_iteration(adata, batch_key, algo_name, seed):
    """Iterates through L1 clusters for sub-clustering without re-loading data."""
    l2_results = []
    l1_clusters = adata.obs[f'leiden_{algo_name}'].unique()

    for cid in l1_clusters:
        logging.info(f"[L2-{algo_name}] Processing Cluster {cid}")
        mask = adata.obs[f'leiden_{algo_name}'] == cid
        sub = adata[mask].copy()

        # Recalculate local manifold (essential to see sub-cluster variance)
        snap.pp.select_features(sub, n_features=50000)
        snap.tl.spectral(sub, features='selected', random_state=seed, distance_metric='cosine', weighted_by_sd=True)

        if algo_name == 'harmony':
            snap.pp.harmony(sub, batch=batch_key, use_rep='X_spectral', random_state=seed)
            rep = 'X_spectral_harmony'
        else:
            snap.pp.mnc_correct(sub, batch=batch_key, use_rep='X_spectral')
            rep = 'X_spectral_mnc'

        snap.pp.knn(sub, n_neighbors=30, use_rep=rep, random_state=seed)
        snap.tl.umap(sub, use_rep=rep, random_state=seed)
        snap.tl.leiden(sub, resolution=0.5, key_added='L2_cluster', random_state=seed)

        # Prepare subset of obs for joining
        sub_obs = sub.obs[['L2_cluster']].copy()
        sub_obs[f'L2_{algo_name}_cluster'] = f"{cid}_" + sub_obs['L2_cluster'].astype(str)
        sub_obs[f'L2_{algo_name}_umap_1'] = sub.obsm['X_umap'][:, 0]
        sub_obs[f'L2_{algo_name}_umap_2'] = sub.obsm['X_umap'][:, 1]

        l2_results.append(sub_obs[[f'L2_{algo_name}_cluster', f'L2_{algo_name}_umap_1', f'L2_{algo_name}_umap_2']])
        del sub

    return pd.concat(l2_results)