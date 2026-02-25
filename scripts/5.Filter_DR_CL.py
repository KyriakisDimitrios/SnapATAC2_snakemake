import numpy as np
import pandas as pd
import sys
import logging
import time
import snapatac2 as snap
from dr_utils import load_data_to_memory
from clustering_utils import run_l2_iteration

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    stream=sys.stderr
)

# --- 0. Setup & Args ---
random_state = 20120224
start_time = time.time()

try:
    args = sys.argv[1:]
    (h5ad_input, batch_var, n_features,max_iter, n_iter,exclude_clusters_str, png_eigen, png_pre,
     png_mnc_s, png_mnc_l, png_harm_s, png_harm_l,
     flag_output, h5ad_output) = args

    n_features, max_iter, n_iter = int(n_features), int(max_iter), int(n_iter)
except Exception as e:
    logging.error(f"Argument Error: {e}")
    sys.exit(1)

exclude_list = [x.strip() for x in exclude_clusters_str.split(',')]
logging.info(f"Filtering OUT leiden_mnc clusters: {exclude_list}")

# --- 1. Load & Global Spectral ---
dat_mem = load_data_to_memory(h5ad_input)
mask = dat_mem.obs['leiden_mnc'].isin(exclude_list)
adata = dat_mem[~mask].copy()
del dat_mem
snap.pp.select_features(adata, n_features=n_features)
snap.tl.spectral(adata, features='selected',
                 random_state=random_state,
                 distance_metric='cosine',
                 weighted_by_sd=True,
                 inplace=True)
snap.pl.spectral_eigenvalues(adata, out_file=png_eigen)

# --- 2. Pre-Correction UMAP ---
snap.tl.umap(adata, use_rep='X_spectral', key_added='umap', random_state=random_state)
snap.pl.umap(adata, color="sample", use_rep='X_umap', out_file=png_pre, height=500, width=500, show=False)

# --- 3. Harmony Block (Global) ---
logging.info("Running Harmony batch correction...")
snap.pp.harmony(adata, batch=batch_var,
                use_rep='X_spectral',
                key_added="X_spectral_harmony",
                max_iter_harmony=max_iter,
                random_state=random_state)
logging.info("Running UMAP on Harmony...")
snap.tl.umap(adata, use_rep="X_spectral_harmony", key_added="umap_harmony", random_state=random_state)
snap.pp.knn(adata, n_neighbors=50, use_rep='X_spectral_harmony', random_state=random_state)
snap.tl.leiden(adata, key_added='leiden_harmony', resolution=2)

snap.pl.umap(adata, color="sample", use_rep="X_umap_harmony", out_file=png_harm_s, height=500, width=500, show=False)
snap.pl.umap(adata, color='leiden_harmony', use_rep="X_umap_harmony", out_file=png_harm_l, height=500, width=500,
             show=False)

# --- 4. MNC Block (Global) ---
logging.info("Running MNC...")
snap.pp.mnc_correct(adata, batch=batch_var, n_iter=n_iter, use_rep='X_spectral', key_added="X_spectral_mnc")
snap.tl.umap(adata, use_rep="X_spectral_mnc", key_added="umap_mnc", random_state=random_state)
snap.pp.knn(adata, n_neighbors=50, use_rep='X_spectral_mnc', random_state=random_state)
snap.tl.leiden(adata, key_added='leiden_mnc', resolution=2)

snap.pl.umap(adata, color="sample", use_rep="X_umap_mnc", out_file=png_mnc_s, height=500, width=500, show=False)
snap.pl.umap(adata, color='leiden_mnc', use_rep="X_umap_mnc", out_file=png_mnc_l, height=500, width=500, show=False)

# --- 5. Iterative Sub-Clustering ---
logging.info("Starting L2 Iterations...")
h_l2 = run_l2_iteration(adata, batch_var, 'harmony', random_state)
m_l2 = run_l2_iteration(adata, batch_var, 'mnc', random_state)

# --- 6. Save & Finish ---
# Merge L2 results back into adata.obs
logging.info("Integrating L2 results into AnnData...")
adata.obs = adata.obs.merge(h_l2, left_index=True, right_index=True, how='left')
adata.obs = adata.obs.merge(m_l2, left_index=True, right_index=True, how='left')

# # Combine all iterative metadata into one CSV
# final_meta = adata.obs.join(h_l2).join(m_l2)
# final_meta.to_csv(csv_iterative)

logging.info(f"Saving AnnData to {h5ad_output}...")
adata.write(h5ad_output, compression="gzip")

with open(flag_output, 'w') as f:
    f.write("done\n")

logging.info(f"Completed Successfully in {time.time() - start_time:.2f}s")