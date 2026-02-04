import sys
import os
import logging
import time
import pandas as pd
import snapatac2 as snap
import scanpy as sc

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

start_time = time.time()
logging.info("Started: Constrain Data (Strict Whitelist)")

try:
    input_file = sys.argv[1]
    sample_name = sys.argv[2]
    annotation_gff3_file = sys.argv[3]
    bin_size = int(sys.argv[4])
    output_file = sys.argv[5]
except IndexError:
    logging.error("Insufficient arguments.")
    sys.exit(1)

logging.info(f"Processing Sample: {sample_name}")


# --- 1. Load Raw Data ---
adata = snap.read(input_file, backed='r', backend='hdf5')
adata_copy = adata.copy(filename=output_file, backend=None)
adata.close()
print(f"Original cells: {adata_copy.n_obs}")


# --- 2. Compute TSS enrichment ---
snap.metrics.tsse(adata_copy, annotation_gff3_file)
logging.info('Completed: Compute TSS enrichment')


# --- 3. add_tile_matrix ---
snap.pp.add_tile_matrix(adata_copy, bin_size=bin_size)
logging.info('Completed: Add Tile Matrix')

# --- 4. Subset & Transfer Metadata ---
# Subset the AnnData
mask = adata_copy.obs['class'] != 'nan'
cells_to_keep = pd.Index(adata_copy.obs_names)[mask]
adata_copy.subset(obs_indices=cells_to_keep, inplace=True)
# Verify
print(f"Remaining cells: {adata_copy.n_obs}")


adata_copy.close()
end_time = time.time()
logging.info(f"Completed in {end_time - start_time:.2f}s")