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
logging.info("Started: Add Metadata")

# --- Inputs from Snakemake ---
try:
    # Unpack all arguments: input files first, then the last two are outputs
    h5ad_input, metadata_path,h5ad_output, flag_file = sys.argv[1:]
except ValueError:
    logging.error("Not enough arguments provided.")
    sys.exit(1)

logging.info(f"H5ad input destination: {h5ad_input}")
logging.info(f"Metadata file: {metadata_path}")
logging.info(f"H5ad output file: {h5ad_output}")
logging.info(f"Flag file: {flag_file}")


# Import Metadata
meta_df = pd.read_csv(metadata_path, index_col='Unnamed: 0')
meta_df['CellID'] = meta_df.index
meta_subset = meta_df.reindex(dat.obs_names)

# Create AnnDataSet
try:
    adata = snap.read(h5ad_input, backed="r", backend=None)
    logging.info(f'Number of cells: {adata.n_obs}')
    logging.info(f'Number of unique barcodes: {np.unique(adata.obs_names).size}')

    adata.close()



except Exception as e:
    logging.error(f"Merge failed: {e}")
    sys.exit(1)

end_time = time.time()
duration = end_time - start_time
logging.info(f"Completed: Merge AnnData")
logging.info(f"Total running time: {duration:.2f} seconds")