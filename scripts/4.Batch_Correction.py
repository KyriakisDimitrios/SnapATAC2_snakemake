import numpy as np
import pandas as pd
import sys
import os
import snapatac2 as snap
import multiprocessing as mp
import logging
import time

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

start_time = time.time()
logging.info("Started: Batch_Correction")

# Inputs from Snakemake
try:
    (h5ad_input,
     batch_var,
     path_to_blacklist,
     n_features,
     max_iter,
     n_iter,
     png_eigenvalue,
     png_sample,
     png_aft_sample,
     png_aft_leiden,
     flag_output,
     h5ad_output) = sys.argv[1:]

    n_features = int(n_features)
    max_iter = int(max_iter)
    n_iter = int(n_iter)
except (IndexError, ValueError):
    logging.error("Not enough arguments provided or incorrect argument types.")
    sys.exit(1)

logging.info(f"Dataset: {h5ad_input}")
logging.info(f"Batch to correct: {batch_var}")
logging.info(f"Blacklist: {path_to_blacklist}")
logging.info(f"Number of features: {n_features}")
logging.info(f"Max iterations: {max_iter}")
logging.info(f"Number of iterations: {n_iter}")
logging.info(f"Eigenvalue PNG: {png_eigenvalue}")
logging.info(f"Sample PNG: {png_sample}")
logging.info(f"Sample AFT PNG: {png_aft_sample}")
logging.info(f"Leiden AFT PNG: {png_aft_leiden}")
logging.info(f"Flag output: {flag_output}")
logging.info(f"H5AD output: {h5ad_output}")



dat = snap.read_dataset(h5ad_input)
# 2. Convert to a single AnnData object (Loads into RAM)
adata = dat.to_adata()
dat.close()
# 3. Write to a single .h5ad file
adata.write(h5ad_output)

dat = snap.read_dataset(h5ad_output)

# select features
logging.info("Selecting features...")
snap.pp.select_features(dat, n_features=n_features,
                        blacklist=path_to_blacklist,
                        max_iter=max_iter,
                        inplace=True,
                        n_jobs=mp.cpu_count() - 2)

logging.info("Running spectral embedding...")
snap.tl.spectral(dat,
                 features='selected',
                 random_state=0,
                 sample_size=None,
                 sample_method='random',
                 chunk_size=20000,
                 distance_metric='cosine',
                 weighted_by_sd=True,
                 feature_weights=None,
                 inplace=True)

logging.info("Running UMAP (pre-correction)...")
snap.tl.umap(dat, n_comps=2, use_dims=None, use_rep='X_spectral', key_added='umap',
             random_state=0, inplace=True)

snap.pl.umap(dat,
             color="sample",
             use_rep='X_umap',
             out_file=png_sample)

# batch correction
logging.info("Running MNC batch correction by sample...")
snap.pp.mnc_correct(dat,
                    batch=batch_var,
                    n_neighbors=5,
                    n_clusters=40,
                    n_iter=n_iter,
                    use_rep='X_spectral',
                    use_dims=None,
                    key_added="X_spectral_mnc_sample_region",
                    inplace=True)

logging.info("Running UMAP (post-correction)...")
snap.tl.umap(dat, n_comps=2, use_dims=None,
             use_rep="X_spectral_mnc_sample_region",
             key_added="umap_mnc_sample_region",
             random_state=0, inplace=True)

snap.pl.umap(dat,
             color="sample",
             use_rep="X_umap_mnc_sample_region",
             out_file=png_aft_sample)

logging.info("Computing KNN graph and Leiden clustering...")
snap.pp.knn(dat,
            n_neighbors=50,
            use_dims=None,
            use_rep='X_spectral_mnc_sample_region',
            method='kdtree',
            inplace=True,
            random_state=0)

snap.tl.leiden(dat)

logging.info("Saving post-correction Leiden UMAP...")
snap.pl.umap(dat, color='leiden', use_rep="X_umap_mnc_sample_region", out_file=png_aft_leiden)

dat.close()

with open(flag_output, 'w') as f:
    pass

end_time = time.time()
duration = end_time - start_time
logging.info("Finished: Batch Correction")
logging.info(f"Total running time: {duration:.2f} seconds")