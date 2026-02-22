import sys
import os
import time
import logging
import functools
import numpy as np
import pandas as pd
import snapatac2 as snap

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    stream=sys.stderr
)

start_time = time.time()
logging.info("Started: Merge AnnData & Consensus Union")

# --- Argument Parsing ---
try:
    annotation_gff3_file = sys.argv[1]
    AnnDataSet_path = sys.argv[2]
    flag_file = sys.argv[3]
    outliers_csv = sys.argv[4]
    n_samples = int(sys.argv[5])

    # Slice the sys.argv list using the known number of samples
    processed_adatas = sys.argv[6: 6 + n_samples]
    features_txt = sys.argv[6 + n_samples: 6 + 2 * n_samples]

except Exception as e:
    logging.error(f"Argument parsing failed: {e}")
    sys.exit(1)

logging.info(f"Initial samples provided: {n_samples}")

# --- Outlier Exclusion (The "Ghost Merge" Bypass) ---
# Assuming the outliers file is in standard results dir. Adjust if necessary.
# We infer the directory from AnnDataSet_path (e.g., results/standard/merge/...)
# base_results_dir = os.path.dirname(os.path.dirname(AnnDataSet_path))
# outliers_csv = os.path.join(base_results_dir, "project_possible_outliers.csv")

bad_samples = set()
if os.path.exists(outliers_csv):
    try:
        bad_df = pd.read_csv(outliers_csv)
        bad_samples = set(bad_df['sample'].astype(str).tolist())
        logging.info(f"Loaded {len(bad_samples)} fatal outliers to exclude from merge.")
    except Exception as e:
        logging.warning(f"Could not load outliers CSV. Error: {e}")

sample_names = [os.path.splitext(os.path.basename(path))[0] for path in processed_adatas]

valid_names = []
valid_adatas = []
valid_features = []

# Filter both the adatas and the feature text files simultaneously
for name, adata_path, feat_path in zip(sample_names, processed_adatas, features_txt):
    if name in bad_samples:
        logging.warning(f"SKIPPING: {name} (Flagged as Outlier)")
        continue
    valid_names.append(name)
    valid_adatas.append(adata_path)
    valid_features.append(feat_path)

logging.info(f"Successfully retained {len(valid_names)} healthy samples for processing.")

# --- 1. Compute Consensus Feature Union ---
logging.info("Calculating consensus feature union across healthy samples...")
try:
    # Read each text file directly into a Python set for O(1) lookups
    all_feature_sets = [set(np.loadtxt(f, dtype=str)) for f in valid_features]

    # Perform rapid mathematical union
    master_features = functools.reduce(set.union, all_feature_sets)
    master_features_list = list(master_features)
    logging.info(f"Consensus Union complete: {len(master_features_list)} unique features identified.")
except Exception as e:
    logging.error(f"Feature union failed: {e}")
    sys.exit(1)

# --- 2. Build AnnDataSet ---
filtered_adatas = []
for name, file in zip(valid_names, valid_adatas):
    adata = snap.read(file, backed="r", backend=None)
    filtered_adatas.append((name, adata))

try:
    # This creates the virtual concatenation
    dat = snap.AnnDataSet(
        adatas=filtered_adatas,
        filename=AnnDataSet_path
    )

    # Generate boolean mask using O(1) set lookups
    # (using the 'master_features' set we calculated above)
    mask = [feat in master_features for feat in dat.var_names]

    # Assign the native python list back to the dataframe
    dat.var['selected'] = mask

    logging.info("Applied consensus feature mask to AnnDataSet.")

    # Add TSS/etc enrichment
    dat.obs['n_fragment'] = dat.adatas.obs['n_fragment']
    dat.obs['frac_dup'] = dat.adatas.obs['frac_dup']
    dat.obs['frac_mito'] = dat.adatas.obs['frac_mito']
    dat.obs['tsse'] = dat.adatas.obs['tsse']
    dat.obs['blacklist_frac'] = dat.adatas.obs['blacklist_frac']

    logging.info('Completed: Transferred observation metadata.')
    logging.info(f'Number of cells: {dat.n_obs}')
    logging.info(f'Number of unique barcodes: {np.unique(dat.obs_names).size}')

    assert dat.n_obs == np.unique(dat.obs_names).size

    dat.close()

    # Create standard execution flag
    with open(flag_file, 'w') as f:
        f.write("done\n")

except Exception as e:
    logging.error(f"Merge failed: {e}")
    sys.exit(1)

duration = time.time() - start_time
logging.info(f"Completed: Merge AnnData in {duration:.2f} seconds")

# import sys
# import os
# import time
# import logging
# import functools
# import numpy as np
# import snapatac2 as snap
#
# # --- Logger Configuration ---
# logging.basicConfig(
#     format='%(asctime)s | %(levelname)-7s | %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     level=logging.INFO,
#     stream=sys.stderr
# )
#
# start_time = time.time()
# logging.info("Started: Merge AnnData & Consensus Union")
#
# # --- Argument Parsing ---
# try:
#     annotation_gff3_file = sys.argv[1]
#     AnnDataSet_path = sys.argv[2]
#     flag_file = sys.argv[3]
#     n_samples = int(sys.argv[4])
#
#     # Slice the sys.argv list using the known number of samples
#     processed_adatas = sys.argv[5: 5 + n_samples]
#     features_txt = sys.argv[5 + n_samples: 5 + 2 * n_samples]
#
# except Exception as e:
#     logging.error(f"Argument parsing failed: {e}")
#     sys.exit(1)
#
# logging.info(f"Merging {n_samples} samples...")
#
# # --- 1. Compute Consensus Feature Union ---
# logging.info("Calculating consensus feature union across all samples...")
# try:
#     # Read each text file directly into a Python set for O(1) lookups
#     all_feature_sets = [set(np.loadtxt(f, dtype=str)) for f in features_txt]
#
#     # Perform rapid mathematical union
#     master_features = functools.reduce(set.union, all_feature_sets)
#     master_features_list = list(master_features)
#     logging.info(f"Consensus Union complete: {len(master_features_list)} unique features identified.")
# except Exception as e:
#     logging.error(f"Feature union failed: {e}")
#     sys.exit(1)
#
# # --- 2. Build AnnDataSet ---
# sample_names = [os.path.splitext(os.path.basename(path))[0] for path in processed_adatas]
# filtered_adatas = []
#
# for name, file in zip(sample_names, processed_adatas):
#     adata = snap.read(file, backed="r", backend=None)
#     filtered_adatas.append((name, adata))
#
# try:
#     # This creates the virtual concatenation
#     dat = snap.AnnDataSet(
#         adatas=filtered_adatas,
#         filename=AnnDataSet_path
#     )
#
#     # Apply the master feature mask directly to the merged object's var dataframe
#     # dat.var['selected'] = dat.var_names.isin(master_features_list)
#
#     # Generate boolean mask using O(1) set lookups
#     # (using the 'master_features' set we calculated above)
#     mask = [feat in master_features for feat in dat.var_names]
#     # Assign the native python list back to the dataframe
#     dat.var['selected'] = mask
#
#     logging.info("Applied consensus feature mask to AnnDataSet.")
#
#     # Add TSS/etc enrichment
#     dat.obs['n_fragment'] = dat.adatas.obs['n_fragment']
#     dat.obs['frac_dup'] = dat.adatas.obs['frac_dup']
#     dat.obs['frac_mito'] = dat.adatas.obs['frac_mito']
#     dat.obs['tsse'] = dat.adatas.obs['tsse']
#     dat.obs['blacklist_frac'] = dat.adatas.obs['blacklist_frac']
#
#     logging.info('Completed: Transferred observation metadata.')
#     logging.info(f'Number of cells: {dat.n_obs}')
#     logging.info(f'Number of unique barcodes: {np.unique(dat.obs_names).size}')
#
#     assert dat.n_obs == np.unique(dat.obs_names).size
#
#     dat.close()
#
#     # Create standard execution flag
#     with open(flag_file, 'w') as f:
#         f.write("done\n")
#
# except Exception as e:
#     logging.error(f"Merge failed: {e}")
#     sys.exit(1)
#
# duration = time.time() - start_time
# logging.info(f"Completed: Merge AnnData in {duration:.2f} seconds")
#


# import numpy as np
# import pandas as pd
# import sys
# import os
# import snapatac2 as snap
# import logging
# import time
#
# # --- Logger Configuration ---
# logging.basicConfig(
#     format='%(asctime)s | %(levelname)-7s | %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     level=logging.INFO,
#     handlers=[logging.StreamHandler(sys.stderr)]
# )
#
# start_time = time.time()
# logging.info("Started: Merge AnnData")
#
# # --- Inputs from Snakemake ---
# # --- Inputs from Snakemake ---
# try:
#     # Unpack all arguments
#     annotation_gff3_file = sys.argv[1]
#     AnnDataSet_path = sys.argv[2]      # <--- FIXED: Removed the trailing comma
#     flag_file = sys.argv[3]
#     processed_adatas = sys.argv[4:]    # Captures all remaining file paths
# except Exception as e:
#     # Changed to catch generic Exception to see the real error if it happens again
#     logging.error(f"Argument parsing failed: {e}")
#     sys.exit(1)
#
#
# logging.info(f"Number of input files: {len(processed_adatas)}")
# logging.info(f"Annotation gff3 file: {annotation_gff3_file}")
# logging.info(f"AnnDataSet destination: {AnnDataSet_path}")
# logging.info(f"Flag file: {flag_file}")
#
# # Extract sample names from filenames
# sample_names = [os.path.splitext(os.path.basename(path))[0] for path in processed_adatas]
# logging.info(f"Sample names: {sample_names}")
#
# filtered_adatas = []
#
# # Loop through samples
# for name, file in zip(sample_names, processed_adatas):
#     # Read in backed mode
#     adata = snap.read(file, backed="r", backend=None)
#
#     # Original logic: (name, adata)
#     filtered_adatas.append((name, adata))
#
# # Create AnnDataSet
# try:
#     dat = snap.AnnDataSet(
#         adatas=filtered_adatas,
#         filename=AnnDataSet_path
#     )
#
#
#     # Add TSS/etc enrichment
#     dat.obs['n_fragment'] = dat.adatas.obs['n_fragment']
#     dat.obs['frac_dup'] = dat.adatas.obs['frac_dup']
#     dat.obs['frac_mito'] = dat.adatas.obs['frac_mito']
#     dat.obs['tsse'] = dat.adatas.obs['tsse']
#     logging.info('Completed: Add TSS enrichment/n_fragment')
#
#     logging.info(f'Number of cells: {dat.n_obs}')
#     logging.info(f'Number of unique barcodes: {np.unique(dat.obs_names).size}')
#
#
#     # Validation
#     assert dat.n_obs == np.unique(dat.obs_names).size
#
#     dat.close()
#
#     # Create flag file
#     with open(flag_file, 'w') as f:
#         f.write("done")
#
# except Exception as e:
#     logging.error(f"Merge failed: {e}")
#     sys.exit(1)
#
# end_time = time.time()
# duration = end_time - start_time
# logging.info(f"Completed: Merge AnnData")
# logging.info(f"Total running time: {duration:.2f} seconds")