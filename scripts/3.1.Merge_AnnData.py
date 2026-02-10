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
# --- Inputs from Snakemake ---
try:
    # Unpack all arguments
    annotation_gff3_file = sys.argv[1]
    AnnDataSet_path = sys.argv[2]      # <--- FIXED: Removed the trailing comma
    flag_file = sys.argv[3]
    processed_adatas = sys.argv[4:]    # Captures all remaining file paths
except Exception as e:
    # Changed to catch generic Exception to see the real error if it happens again
    logging.error(f"Argument parsing failed: {e}")
    sys.exit(1)


logging.info(f"Number of input files: {len(processed_adatas)}")
logging.info(f"Annotation gff3 file: {annotation_gff3_file}")
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


    # Add TSS/etc enrichment
    dat.obs['n_fragment'] = dat.adatas.obs['n_fragment']
    dat.obs['frac_dup'] = dat.adatas.obs['frac_dup']
    dat.obs['frac_mito'] = dat.adatas.obs['frac_mito']
    dat.obs['tsse'] = dat.adatas.obs['tsse']
    logging.info('Completed: Add TSS enrichment/n_fragment')

    logging.info(f'Number of cells: {dat.n_obs}')
    logging.info(f'Number of unique barcodes: {np.unique(dat.obs_names).size}')


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