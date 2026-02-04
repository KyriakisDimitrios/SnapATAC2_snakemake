import sys
import os
import logging
import time
import pandas as pd
import numpy as np
import snapatac2 as snap

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

start_time = time.time()
logging.info("Started: Constrain Data (Dynamic Grouping)")

# --- ARGUMENT PARSING ---
try:
    input_file = sys.argv[1]
    sample_name = sys.argv[2]
    annotation_gff3_file = sys.argv[3]
    # 4. Target Subtypes (Comma-separated string)
    target_subtypes_str = sys.argv[4]
    target_subtypes = target_subtypes_str.split(',')
    # 5. Group By Column
    group_by_col = sys.argv[5]
    bin_size = int(sys.argv[6])
    output_file = sys.argv[7]
except IndexError:
    logging.error("Insufficient arguments. Expected: input sample annotation targets group_col bin_size output")
    sys.exit(1)

logging.info(f"Processing Sample: {sample_name}")
logging.info(f"Filtering on column: '{group_by_col}'")
logging.info(f"Targets: {target_subtypes}")

# --- 1. Load Raw Data ---
adata = snap.read(input_file, backed='r', backend='hdf5')
adata_copy = adata.copy(filename=output_file, backend=None)
adata.close()

# --- 2. Compute TSS enrichment ---
snap.metrics.tsse(adata_copy, annotation_gff3_file)
logging.info('Completed: Compute TSS enrichment')

# --- 3. Add Tile Matrix ---
snap.pp.add_tile_matrix(adata_copy, bin_size=bin_size)
logging.info('Completed: Add Tile Matrix')

# --- 4. Subset Cell Types ---
# 1. Get the column data into memory as a list/array
subtypes = list(adata_copy.obs[group_by_col])
# 2. Create the mask using NumPy (faster/safer than pandas for raw lists)
mask = np.isin(subtypes, target_subtypes)
# 3. Get the specific Cell IDs to keep
cells_to_keep = pd.Index(adata_copy.obs_names)[mask]

print(f"Original: {adata_copy.n_obs}")
print(f"Filtered: {len(cells_to_keep)}")

# 5. Apply the subset
adata_copy.subset(obs_indices=cells_to_keep, inplace=True)

if len(cells_to_keep) == 0:
    logging.error(f"No cells found matching {target_subtypes} in {group_by_col}.")
    sys.exit(1)

# 4. Apply Subset
adata_copy.subset(obs_indices=cells_to_keep, inplace=True)

# Save and Close
adata_copy.close()
end_time = time.time()
logging.info(f"Completed in {end_time - start_time:.2f}s")