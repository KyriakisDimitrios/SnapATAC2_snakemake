import sys
import os
import logging
import pandas as pd
import scanpy as sc

# --- 1. Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(message)s')

# Get arguments from command line
input_h5ad  = sys.argv[1]
rush_meta_path = sys.argv[2]
va_meta_path   = sys.argv[3]
output_h5ad = sys.argv[4]

logging.info(f"Processing: {os.path.basename(input_h5ad)}")

# --- 2. Load the Data ---
logging.info("Loading AnnData...")
adata = sc.read_h5ad(input_h5ad)

# --- 3. Prepare the Metadata (The Lookup Table) ---
logging.info("Loading Metadata...")

# A. Load Rush Metadata
df_rush = pd.read_csv(rush_meta_path, sep="\t")
# Create 'TestID' = "Rush_054" (Clean up spaces and names)
df_rush['clean_name'] = df_rush['sample_name'].astype(str).str.replace(r'^[^_]+_', '', regex=True)
df_rush['TestID'] = df_rush['BB'].astype(str).str.strip() + "_" + df_rush['clean_name']
df_rush['TestID'] = df_rush['TestID'].str.replace(' ', '') # Remove spaces

# B. Load VA Metadata
df_va = pd.read_csv(va_meta_path, sep="\t")
# Create 'TestID' = "VA-2302" (Standardize format)
df_va['TestID'] = df_va['sample_name'].astype(str).str.replace('_', '-', regex=False)

# C. Combine them
lookup_table = pd.concat([df_va, df_rush], ignore_index=True)
# Remove duplicates just in case
lookup_table = lookup_table.drop_duplicates(subset=['TestID'])
# Set the index so we can look things up easily
lookup_table = lookup_table.set_index('TestID')

logging.info(f"Metadata loaded. Found {len(lookup_table)} samples in the table.")

# --- 4. Map Data to Cells ---
logging.info("Mapping metadata to cells...")

# We create a new column 'clean_sample_id' in adata to match the format of our table
# e.g., "Rush-54_multi" -> "Rush_054"

# Step A: Get raw sample names from the adata
# We make a copy to work on
adata.obs['donor'] = adata.obs['sample'].astype(str)
adata.obs['clean_sample_id'] = adata.obs['sample'].astype(str)

# Step B: Fix formatting to match the Lookup Table
# 1. Remove '_multi'
adata.obs['clean_sample_id'] = adata.obs['clean_sample_id'].str.replace(r'_multi$', '', regex=True)
# 2. Fix Rush IDs (Rush-54 -> Rush_054)
#    We split by '-', take the number, pad it with zeros, and put it back together
#    Note: This specific line handles the "Rush-54" to "Rush_054" conversion
def fix_rush_id(x):
    if "Rush" in x and "-" in x:
        parts = x.split("-")
        if len(parts) == 2 and parts[1].isdigit():
            return f"Rush_{int(parts[1]):03d}"
    return x.replace("-", "_") if "Rush" in x else x

adata.obs['clean_sample_id'] = adata.obs['clean_sample_id'].apply(fix_rush_id)

# Step C: The Merge (This is the magic part)
# We find the intersection of sample IDs
valid_ids = adata.obs['clean_sample_id'].unique()
logging.info(f"Found {len(valid_ids)} unique samples in the data.")

# Create a dictionary for every column in the metadata
# e.g. map_dict = {'Rush_054': 'Male', 'VA-2302': 'Female'}
for col in lookup_table.columns:
    # Create the map
    map_dict = lookup_table[col].to_dict()
    # Apply it to the cells
    adata.obs[col] = adata.obs['clean_sample_id'].map(map_dict)

# --- 5. Save ---
logging.info("Sanitizing metadata columns...")
for col in adata.obs.columns:
    # Force convert all object columns to string to handle NaNs/Integers mixed in
    if adata.obs[col].dtype == 'object':
        adata.obs[col] = adata.obs[col].astype(str)
logging.info(f"Saving to {output_h5ad}...")

# Ensure the folder exists
output_dir = os.path.dirname(output_h5ad)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)

adata.write(output_h5ad, compression="gzip")
logging.info("Done.")