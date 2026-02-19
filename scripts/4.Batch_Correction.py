import numpy as np
import pandas as pd
import sys
import os
import snapatac2 as snap
import multiprocessing as mp
import logging
import time
random_state = 20120224

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

start_time = time.time()
logging.info("Started: Batch_Correction")

# --- Argument Parsing ---
# We expect exactly 14 arguments from the Snakemake rule
try:
    args = sys.argv[1:]
    if len(args) != 14:
        raise ValueError(f"Expected 14 arguments, got {len(args)}")

    (h5ad_input,
     batch_var,
     path_to_blacklist,
     n_features,
     max_iter,
     n_iter,
     png_eigen,       # 1
     png_pre,         # 2
     png_mnc_s,       # 3
     png_mnc_l,       # 4
     png_harm_s,      # 5
     png_harm_l,      # 6
     flag_output,
     h5ad_output) = args

    n_features = int(n_features)
    max_iter = int(max_iter)
    n_iter = int(n_iter)

except Exception as e:
    logging.error(f"Argument Error: {e}")
    sys.exit(1)

logging.info(f"Dataset: {h5ad_input}")
logging.info(f"Batch Variable: {batch_var}")

# --- 1. Load Data ---
# Load as a backed dataset to save memory
dat = snap.read_dataset(h5ad_input)
# dat = snap.read(h5ad_input, backed=None)
# --- 2. Feature Selection ---
logging.info("Selecting features...")
snap.pp.select_features(dat, n_features=n_features,
                        blacklist=path_to_blacklist,
                        max_iter=max_iter,
                        inplace=True,
                        n_jobs=mp.cpu_count() - 2)

# --- 3. Spectral Embedding ---
logging.info("Running spectral embedding...")
snap.tl.spectral(dat,
                 features='selected',
                 random_state=random_state,
                 distance_metric='cosine',
                 weighted_by_sd=True,
                 inplace=True)

# Plot 1: Eigenvalues
snap.pl.spectral_eigenvalues(dat, out_file=png_eigen)

# --- 4. Pre-Correction UMAP ---
logging.info("Running UMAP (Pre-correction)...")
snap.tl.umap(dat, n_comps=2, use_rep='X_spectral', key_added='umap', random_state=0, inplace=True)

# Plot 2: Pre-Correction Sample
snap.pl.umap(dat, color="sample", use_rep='X_umap', out_file=png_pre, height=500, width=500, show=False)


# --- 5. Harmony Correction ---
logging.info("Running Harmony batch correction...")
snap.pp.harmony(dat, batch=batch_var,
                  use_rep='X_spectral',
                  key_added="X_spectral_harmony",
                  max_iter_harmony=20,
                  random_state=random_state)
logging.info("Running UMAP on Harmony...")
snap.tl.umap(dat, n_comps=2,
                 use_rep="X_spectral_harmony",
                 key_added="umap_harmony",
                 random_state=random_state,
                 inplace=True)
logging.info("Clustering Harmony...")
snap.pp.knn(dat, n_neighbors=50,
            use_rep='X_spectral_harmony',
            method='kdtree', inplace=True, random_state=random_state)
snap.tl.leiden(dat, key_added='leiden_harmony',resolution=2)


# Plot 5 & 6: Harmony Sample + Harmony Leiden
snap.pl.umap(dat, color="sample", use_rep="X_umap_harmony", out_file=png_harm_s, height=500, width=500, show=False)
snap.pl.umap(dat, color='leiden_harmony', use_rep="X_umap_harmony", out_file=png_harm_l, height=500, width=500, show=False)


# --- 5. MNC Correction ---
logging.info("Running MNC batch correction...")
snap.pp.mnc_correct(dat,
                    batch=batch_var,
                    n_neighbors=5,
                    n_clusters=40,
                    n_iter=n_iter,
                    use_rep='X_spectral',
                    key_added="X_spectral_mnc",
                    inplace=True)

logging.info("Running UMAP on MNC...")
snap.tl.umap(dat, n_comps=2,
             use_rep="X_spectral_mnc",
             random_state=random_state,
             key_added="umap_mnc",
             inplace=True)

# Cluster MNC
logging.info("Clustering MNC...")
snap.pp.knn(dat, n_neighbors=50, use_rep='X_spectral_mnc', method='kdtree', inplace=True, random_state=0)
snap.tl.leiden(dat, key_added='leiden_mnc',resolution=2)

# Plot 3 & 4: MNC Sample + MNC Leiden
snap.pl.umap(dat, color="sample", use_rep="X_umap_mnc", out_file=png_mnc_s, height=500, width=500, show=False)
snap.pl.umap(dat, color='leiden_mnc', use_rep="X_umap_mnc", out_file=png_mnc_l, height=500, width=500, show=False)


# --- 7. Save Final Output ---
logging.info(f"Saving AnnData to {h5ad_output}...")

# Convert to in-memory AnnData to write strictly as .h5ad
# Note: snap.read_dataset returns an AnnDataSet (backed).
# .to_adata() converts it to a standard AnnData object.
out_adata = dat.to_adata()
out_adata.write(h5ad_output, compression="gzip")

dat.close()

# Create Flag
with open(flag_output, 'w') as f:
    f.write("Batch correction completed successfully.")

end_time = time.time()
logging.info(f"Finished. Total time: {end_time - start_time:.2f}s")