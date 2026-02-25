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
    (h5ad_input, batch_var, max_iter, n_iter, png_eigen, png_pre,
     png_mnc_s, png_mnc_l, png_harm_s, png_harm_l,
     flag_output, h5ad_output) = args

    max_iter, n_iter = int(max_iter), int(n_iter)
except Exception as e:
    logging.error(f"Argument Error: {e}")
    sys.exit(1)

# --- 1. Load & Global Spectral ---
dat_mem = load_data_to_memory(h5ad_input)
num_features = dat_mem.var['selected'].sum()
logging.info(f"Using {num_features} pre-calculated stratified consensus features.")

snap.tl.spectral(dat_mem, features='selected',
                 random_state=random_state,
                 distance_metric='cosine',
                 weighted_by_sd=True,
                 inplace=True)
snap.pl.spectral_eigenvalues(dat_mem, out_file=png_eigen)

# --- 2. Pre-Correction UMAP ---
snap.tl.umap(dat_mem, use_rep='X_spectral', key_added='umap', random_state=random_state)
snap.pl.umap(dat_mem, color="sample", use_rep='X_umap', out_file=png_pre, height=500, width=500, show=False)

# --- 3. Harmony Block (Global) ---
logging.info("Running Harmony batch correction...")
snap.pp.harmony(dat_mem, batch=batch_var,
                use_rep='X_spectral',
                key_added="X_spectral_harmony",
                max_iter_harmony=max_iter,
                random_state=random_state)
logging.info("Running UMAP on Harmony...")
snap.tl.umap(dat_mem, use_rep="X_spectral_harmony", key_added="umap_harmony", random_state=random_state)
snap.pp.knn(dat_mem, n_neighbors=50, use_rep='X_spectral_harmony', random_state=random_state)
snap.tl.leiden(dat_mem, key_added='leiden_harmony', resolution=2)

snap.pl.umap(dat_mem, color="sample", use_rep="X_umap_harmony", out_file=png_harm_s, height=500, width=500, show=False)
snap.pl.umap(dat_mem, color='leiden_harmony', use_rep="X_umap_harmony", out_file=png_harm_l, height=500, width=500,
             show=False)

# --- 4. MNC Block (Global) ---
logging.info("Running MNC...")
snap.pp.mnc_correct(dat_mem, batch=batch_var, n_iter=n_iter, use_rep='X_spectral', key_added="X_spectral_mnc")
snap.tl.umap(dat_mem, use_rep="X_spectral_mnc", key_added="umap_mnc", random_state=random_state)
snap.pp.knn(dat_mem, n_neighbors=50, use_rep='X_spectral_mnc', random_state=random_state)
snap.tl.leiden(dat_mem, key_added='leiden_mnc', resolution=2)

snap.pl.umap(dat_mem, color="sample", use_rep="X_umap_mnc", out_file=png_mnc_s, height=500, width=500, show=False)
snap.pl.umap(dat_mem, color='leiden_mnc', use_rep="X_umap_mnc", out_file=png_mnc_l, height=500, width=500, show=False)

# --- 5. Iterative Sub-Clustering ---
logging.info("Starting L2 Iterations...")
h_l2 = run_l2_iteration(dat_mem, batch_var, 'harmony', random_state)
m_l2 = run_l2_iteration(dat_mem, batch_var, 'mnc', random_state)

# --- 6. Save & Finish ---
# Merge L2 results back into dat_mem.obs
logging.info("Integrating L2 results into AnnData...")
dat_mem.obs = dat_mem.obs.merge(h_l2, left_index=True, right_index=True, how='left')
dat_mem.obs = dat_mem.obs.merge(m_l2, left_index=True, right_index=True, how='left')

# # Combine all iterative metadata into one CSV
# final_meta = dat_mem.obs.join(h_l2).join(m_l2)
# final_meta.to_csv(csv_iterative)

logging.info(f"Saving AnnData to {h5ad_output}...")
dat_mem.write(h5ad_output, compression="gzip")

with open(flag_output, 'w') as f:
    f.write("done\n")

logging.info(f"Completed Successfully in {time.time() - start_time:.2f}s")