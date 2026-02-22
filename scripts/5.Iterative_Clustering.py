import os
import sys
import logging
import pandas as pd
import snapatac2 as snap

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)


def _atac_clustering_block(adata, n_features, leiden_res, batch_key, level_prefix, algorithm):
    """Runs ATAC dimensionality reduction with dynamic batch correction."""

    logging.info(f"[{level_prefix}] Selecting top {n_features} features...")
    snap.pp.select_features(adata, n_features=n_features)

    logging.info(f"[{level_prefix}] Running Spectral Embedding...")
    snap.tl.spectral(adata, features='selected')

    # --- Dynamic Batch Correction Switch ---
    if algorithm.lower() == 'harmony':
        logging.info(f"[{level_prefix}] Running Harmony Batch Correction on '{batch_key}'...")
        snap.pp.harmony(adata, batch=batch_key)
        rep_use = 'harmony'
    elif algorithm.lower() == 'mnc':
        logging.info(f"[{level_prefix}] Running MNC Batch Correction on '{batch_key}'...")
        snap.pp.mnc(adata, batch=batch_key)
        rep_use = 'mnc'
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
                              l1_features=100000, l1_res=0.5,
                              l2_features=50000, l2_res=0.8):
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # --- LEVEL 1: Global Clustering ---
    logging.info(f"Loading {input_h5ad} entirely into RAM (backed=None)...")
    adata = snap.read(input_h5ad, backed=None)

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


# --- Argument Parsing ---
if __name__ == "__main__":
    try:
        in_h5ad = sys.argv[1]
        out_csv = sys.argv[2]
        b_key = sys.argv[3]
        algo = sys.argv[4]
    except IndexError:
        logging.error("Usage: script.py <input.h5ad> <output.csv> <batch_key> <algorithm>")
        sys.exit(1)

    iterative_atac_clustering(
        input_h5ad=in_h5ad,
        output_csv=out_csv,
        batch_key=b_key,
        algorithm=algo
    )