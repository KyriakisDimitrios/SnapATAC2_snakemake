import numpy as np
import pandas as pd
import sys
import os
import snapatac2 as snap
import logging
import time
import scanpy as sc
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
    h5ad_input = sys.argv[1]
    rush_meta   = sys.argv[2]
    va_meta = sys.argv[3]
    h5ad_output = sys.argv[4]
except ValueError:
    logging.error("Not enough arguments provided.")
    sys.exit(1)

logging.info(f"H5ad input destination: {h5ad_input}")
logging.info(f"Metadata Rush file: {rush_meta}")
logging.info(f"Metadata VA file: {va_meta}")
logging.info(f"H5ad output file: {h5ad_output}")

# Load Rush
meta_rush = pd.read_csv(rush_meta, sep="\t")
meta_rush['BB'] = meta_rush['BB'].astype(str).str.strip()  # Remove whitespace
cleaned_names = meta_rush['sample_name'].astype(str).str.replace(r'^[^_]+_', '', regex=True)
meta_rush['TestID'] = meta_rush['BB'] + "_" + cleaned_names
# Fix: Ensure Rush IDs are clean (remove spaces like 'Rush _124')
meta_rush['TestID'] = meta_rush['TestID'].str.replace(' ', '')

# Load VA
meta_va = pd.read_csv(va_meta, sep="\t")
meta_va['TestID'] = meta_va['sample_name'].astype(str).str.replace('_', '-', regex=False)

# Combine
combined_meta = pd.concat([meta_va, meta_rush], ignore_index=True)
combined_meta = combined_meta.drop_duplicates(subset=['TestID'])
meta_lookup = combined_meta.set_index('TestID')

logging.info(f"Total metadata rows: {len(meta_lookup)}")
logging.info(f"Sample Check - Rush: 'Rush_054' in metadata? {'Rush_054' in meta_lookup.index}")


# ====================== 3. Logic Function ======================
def get_metadata_row(donor_id, meta_lookup):
    """
    Smart lookup:
    1. Exact Match
    2. Swap Separator (Rush-54 -> Rush_54)
    3. Zero Pad (Rush_54 -> Rush_054)
    """
    # 1. Exact
    if donor_id in meta_lookup.index:
        return meta_lookup.loc[donor_id]

    # 2. Separator Swap
    alt_id = donor_id.replace('-', '_')
    if alt_id in meta_lookup.index:
        return meta_lookup.loc[alt_id]

    # 3. Zero Pad (Specific for Rush)
    if "Rush" in alt_id:
        try:
            parts = alt_id.split('_')
            # If we have [Rush, 54], pad the number
            if len(parts) == 2 and parts[1].isdigit():
                padded_id = f"{parts[0]}_{int(parts[1]):03d}"
                if padded_id in meta_lookup.index:
                    logging.info(f"   (Logic Applied: {donor_id} -> {padded_id})")
                    return meta_lookup.loc[padded_id]
        except:
            pass
    return None


# ====================== 4. Test on File (SnapATAC2) ======================
logging.info(f"\n--- Processing {os.path.basename(h5ad_input)} ---")

# Open in backed mode ('r+') to allow writing metadata without loading matrix
adata = sc.read_h5ad(h5ad_input)

# A. Create Donor Column
if 'sample' in adata.obs.columns:
    # Remove _multi suffix
    adata.obs['donor'] = adata.obs['sample'].astype(str).str.replace(r'_multi$', '', regex=True)
    current_donor = str(adata.obs['donor'].iloc[0])
    logging.info(f"File Donor ID: {current_donor}")

    # B. Find Match
    row = get_metadata_row(current_donor, meta_lookup)

    if row is not None:
        logging.info(f"MATCH FOUND: Linked to {row.name}")

        # C. Write Metadata
        for col in row.index:
            adata.obs[col] = row[col]
        logging.info("Metadata successfully written to adata.obs")
    else:
        logging.info("!! FAILURE: No match found in metadata !!")
else:
    logging.info("Error: 'sample' column missing")

adata.write(h5ad_output)

end_time = time.time()
duration = end_time - start_time
logging.info(f"Completed: Add Metadata")
logging.info(f"Total running time: {duration:.2f} seconds")