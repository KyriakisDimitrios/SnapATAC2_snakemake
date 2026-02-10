import sys
import os
import logging
from unittest.mock import inplace
import pandas as pd
import snapatac2 as snap
import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(message)s')

# --- 0. Parse Arguments ---
input_h5ad = sys.argv[1]
output_h5ad = sys.argv[2]
min_tsse = float(sys.argv[3])
qc_files = sys.argv[4:]

# --- 1. Load QC Stats ---
logging.info("--- Step 1: Loading QC Stats ---")
qc_scores = {}
for fpath in qc_files:
    try:
        fname = os.path.basename(fpath)
        # Clean filename: "Rush-01_raw_qc.txt" -> "Rush-01"
        sample = fname.replace("_raw_qc.txt", "").replace(".raw_qc.txt", "")
        with open(fpath, 'r') as f:
            score = float(f.read().strip().split(',')[0])
            qc_scores[sample] = score
    except:
        pass

# Identify passing samples
valid_samples = {s for s, score in qc_scores.items() if score >= min_tsse}
logging.info(f"Threshold: {min_tsse} | Passing Samples: {len(valid_samples)}")

if not valid_samples:
    logging.error("No samples passed the threshold!")
    sys.exit(1)

# --- 3. Stream Filter ---
logging.info(f"--- Step 2: Filtering from {input_h5ad} ---")

# CRITICAL: Open in 'r+' (Read/Write) mode to allow internal locking during subset
adata = snap.read_dataset(input_h5ad)

# Convert column to list (as requested)
all_samples = list(adata.obs['sample'])

# Clean suffix manually in the list (Remove _multi)

# Create Boolean Mask
# Convert valid_samples to list for numpy
mask = np.isin(all_samples, list(valid_samples))
cells_to_keep = pd.Index(adata.obs_names)[mask]


# WRITE TO NEW FILE
# Using 'out=' creates a fresh, defragmented file (Fixes 'Empty Slot' crash)
logging.info(f"Writing to {output_h5ad}...")
adata.subset(obs_indices=cells_to_keep,out=output_h5ad)

adata.close()
logging.info("Done.")