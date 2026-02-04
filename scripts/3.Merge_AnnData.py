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
logging.info("Started: Merge AnnData")

# --- Inputs from Snakemake ---
try:
    # Unpack all arguments: input files first, then the last two are outputs
    *processed_adatas, metadata_path,AnnDataSet_path, flag_file = sys.argv[1:]
except ValueError:
    logging.error("Not enough arguments provided.")
    sys.exit(1)

logging.info(f"Number of input files: {len(processed_adatas)}")
logging.info(f"AnnDataSet destination: {AnnDataSet_path}")
logging.info(f"Flag file: {flag_file}")

# Extract sample names from filenames
sample_names = [os.path.splitext(os.path.basename(path))[0] for path in processed_adatas]
logging.info(f"Sample names: {sample_names}")

filtered_adatas = []

# Loop through samples
for name, file in zip(sample_names, processed_adatas):
    # Read in backed mode
    adata = snap.read(file, backed="r", backend=None)

    # Original logic: (name, adata)
    filtered_adatas.append((name, adata))

# Create AnnDataSet
try:
    dat = snap.AnnDataSet(
        adatas=filtered_adatas,
        filename=AnnDataSet_path
    )

    logging.info(f'Number of cells: {dat.n_obs}')
    logging.info(f'Number of unique barcodes: {np.unique(dat.obs_names).size}')

    # Import Metadata
    meta_df = pd.read_csv(metadata_path, index_col='Unnamed: 0')
    meta_df['CellID'] = meta_df.index
    meta_subset = meta_df.reindex(dat.obs_names)
    for col in meta_subset.columns:
        series = meta_subset[col]
        if series.dtype == 'object' or series.dtype.name == 'category':
            values = series.astype(str).values
        else:
            values = series.values
        dat.obs[col] = values
    
    # Create unique cell IDs
    unique_cell_ids = [sa + ':' + bc for sa, bc in zip(dat.obs['sample'], dat.obs_names)]
    dat.obs_names = unique_cell_ids

    # Validation
    assert dat.n_obs == np.unique(dat.obs_names).size

    dat.close()

    # Create flag file
    with open(flag_file, 'w') as f:
        f.write("done")

except Exception as e:
    logging.error(f"Merge failed: {e}")
    sys.exit(1)

end_time = time.time()
duration = end_time - start_time
logging.info(f"Completed: Merge AnnData")
logging.info(f"Total running time: {duration:.2f} seconds")