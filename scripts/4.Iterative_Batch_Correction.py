import sys
import os
import warnings

import numpy as np
import pandas as pd

import scanpy as sc
import snapatac2 as snap

import multiprocessing as mp
import logging
import time
import concurrent.futures
from tqdm.auto import tqdm
from dr_utils import load_data_to_memory
from dr_utils import split_to_memory_folds
from dr_utils import mapped_robust_features
from dr_utils import get_hpc_resources

# Hide the specific Polars DeprecationWarning from SnapATAC2
warnings.filterwarnings("ignore", category=DeprecationWarning, message=".*_import_from_c.*")

# --- Logger Configuration ---
# All logs go strictly to stderr as requested
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    stream=sys.stderr
)

random_state = 20120224
start_time = time.time()
logging.info("Started: Batch_Correction")

# --- Argument Parsing ---
# We expect exactly 14 arguments from the Snakemake rule
try:
    args = sys.argv[1:]
    if len(args) != 16:
        raise ValueError(f"Expected 16 arguments, got {len(args)}")

    (h5ad_input,
     batch_var,
     path_to_blacklist,
     n_folds,
     cells_per_sample_fold,
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

    n_folds = int(n_folds)
    cells_per_sample_fold = int(cells_per_sample_fold)
    n_features = int(n_features)
    max_iter = int(max_iter)
    n_iter = int(n_iter)

except Exception as e:
    logging.error(f"Argument Error: {e}")
    sys.exit(1)

logging.info(f"Dataset: {h5ad_input}")
logging.info(f"Batch Variable: {batch_var}")



# --- Execution (Run in its own cell) ---
# 1. Load to RAM
dat_mem = load_data_to_memory(h5ad_input)

# 2. Generate folds
memory_folds = split_to_memory_folds(dat_mem=dat_mem, batch_var = batch_var, n_folds=n_folds, cells_per_sample=cells_per_sample_fold)


# --- Execution Logic ---
n_folds_actual = len(memory_folds)
n_fold_jobs, n_jobs_per_fold = get_hpc_resources(n_folds_actual, n_folds_actual)

robust_features = mapped_robust_features(
    memory_folds,
    n_fold_jobs=n_fold_jobs,
    n_jobs_per_fold=n_jobs_per_fold,
    n_features=n_features,
    sample_size=None # High precision
)

# Apply the mathematically verified features to your master dataset
logging.info("Applying robust features to master dataset...")
dat_mem.var['selected'] = dat_mem.var_names.isin(robust_features)


# --- 3. Spectral Embedding ---
logging.info("Running spectral embedding...")
snap.tl.spectral(dat_mem,
                 features='selected',
                 random_state=random_state,
                 distance_metric='cosine',
                 weighted_by_sd=True,
                 inplace=True)

# Plot 1: Eigenvalues
snap.pl.spectral_eigenvalues(dat_mem, out_file=png_eigen)

# --- 4. Pre-Correction UMAP ---
logging.info("Running UMAP (Pre-correction)...")
snap.tl.umap(dat_mem, n_comps=2, use_rep='X_spectral', key_added='umap', random_state=0, inplace=True)

# Plot 2: Pre-Correction Sample
snap.pl.umap(dat_mem, color="sample", use_rep='X_umap', out_file=png_pre, height=500, width=500, show=False)


# --- 5. Harmony Correction ---
logging.info("Running Harmony batch correction...")
snap.pp.harmony(dat_mem, batch=batch_var,
                  use_rep='X_spectral',
                  key_added="X_spectral_harmony",
                  max_iter_harmony=20,
                  random_state=random_state)
logging.info("Running UMAP on Harmony...")
snap.tl.umap(dat_mem, n_comps=2,
                 use_rep="X_spectral_harmony",
                 key_added="umap_harmony",
                 random_state=random_state,
                 inplace=True)
logging.info("Clustering Harmony...")
snap.pp.knn(dat_mem, n_neighbors=50,
            use_rep='X_spectral_harmony',
            method='kdtree', inplace=True, random_state=random_state)
snap.tl.leiden(dat_mem, key_added='leiden_harmony',resolution=2)


# Plot 5 & 6: Harmony Sample + Harmony Leiden
snap.pl.umap(dat_mem, color="sample", use_rep="X_umap_harmony", out_file=png_harm_s, height=500, width=500, show=False)
snap.pl.umap(dat_mem, color='leiden_harmony', use_rep="X_umap_harmony", out_file=png_harm_l, height=500, width=500, show=False)


dat_mem.write("robust_itr_batch_harmony.h5ad", compression="gzip")


# --- 5. MNC Correction ---
logging.info("Running MNC batch correction...")
snap.pp.mnc_correct(dat_mem,
                    batch=batch_var,
                    n_neighbors=5,
                    n_clusters=40,
                    n_iter=n_iter,
                    use_rep='X_spectral',
                    key_added="X_spectral_mnc",
                    inplace=True)

logging.info("Running UMAP on MNC...")
snap.tl.umap(dat_mem, n_comps=2,
             use_rep="X_spectral_mnc",
             random_state=random_state,
             key_added="umap_mnc",
             inplace=True)

# Cluster MNC
logging.info("Clustering MNC...")
snap.pp.knn(dat_mem, n_neighbors=50, use_rep='X_spectral_mnc', method='kdtree', inplace=True, random_state=0)
snap.tl.leiden(dat_mem, key_added='leiden_mnc',resolution=2)

# Plot 3 & 4: MNC Sample + MNC Leiden
snap.pl.umap(dat_mem, color="sample", use_rep="X_umap_mnc", out_file=png_mnc_s, height=500, width=500, show=False)
snap.pl.umap(dat_mem, color='leiden_mnc', use_rep="X_umap_mnc", out_file=png_mnc_l, height=500, width=500, show=False)


# --- 7. Save Final Output ---
logging.info(f"Saving AnnData to {h5ad_output}...")

logging.info(f"Saving {len(robust_features)} robust features to text file...")
np.savetxt(
    "robust_features_consensus.txt",
    robust_features,
    fmt='%s',
    header="feature_name"

)


dat_mem.write("robust_itr_batch_all.h5ad", compression="gzip")


dat_mem.write(h5ad_output, compression="gzip")

# Create the flag file that Snakemake is looking for
with open(flag_output, 'w') as f:
    f.write("done\n")
logging.info(f"Flag file created at {flag_output}")

logging.info(f"Completed Successful")


end_time = time.time()
logging.info(f"Finished. Total time: {end_time - start_time:.2f}s")