import numpy as np
import pandas as pd
import sys
import os
import snapatac2 as snap
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
logging.info("Started: Clustering")

# --- Argument Parsing ---
try:
    flag_input = sys.argv[1]
    dataset = sys.argv[2]
    resolutions_str = sys.argv[3]
    clustering_mnc = sys.argv[4]
    output_h5ad_path = sys.argv[5]
    flag_output = sys.argv[6]
except IndexError:
    logging.error("Not enough arguments provided.")
    sys.exit(1)

logging.info(f"Dataset: {dataset}")
logging.info(f"Resolutions: {resolutions_str}")
logging.info(f"Output path: {output_h5ad_path}")

# Create output directory
if not os.path.exists(clustering_mnc):
    os.makedirs(clustering_mnc)

# Parse resolutions
resolutions = [float(x.strip()) for x in resolutions_str.split(',')]

# --- Execution ---
dat = snap.read_dataset(dataset)
logging.info("Dataset loaded.")

for resolution in resolutions:
    logging.info(f"Running Leiden clustering at resolution: {resolution}")
    snap.tl.leiden(
        dat,
        resolution=resolution,
        objective_function="modularity",
        min_cluster_size=3,
        n_iterations=-1,
        random_state=2012224,
        key_added=f"leiden_mnc_{resolution}",
        use_leidenalg=False,
        inplace=True
    )

    output_file_path = os.path.join(clustering_mnc, f'Leiden_CL_res_{resolution}.png')
    snap.pl.umap(
        dat,
        color=f"leiden_mnc_{resolution}",
        use_rep='X_umap_mnc_sample_region',
        out_file=output_file_path
    )
    logging.info(f"Saved UMAP for resolution {resolution}")

# Export to AnnData
logging.info("Converting AnnDataSet to AnnData...")

# --- SAFETY WRAPPER START ---
# We wrap ONLY this line so if it fails, the script continues and saves your work.
try:
    dat.obsm['fragment_paired'] = dat.adatas.obsm['fragment_paired']
except Exception as e:
    logging.warning(f"Could not transfer 'fragment_paired'. Skipping to save. Error: {e}")
# --- SAFETY WRAPPER END ---

adat = dat.to_adata()
dat.close()

# Save output
adat.write_h5ad(output_h5ad_path, compression="gzip")
logging.info(f"AnnData saved to {output_h5ad_path}")

# Create flag file
with open(flag_output, 'w') as f:
    pass

end_time = time.time()
duration = end_time - start_time
logging.info("Finished: Clustering")
logging.info(f"Total running time: {duration:.2f} seconds")