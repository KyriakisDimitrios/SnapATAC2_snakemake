import sys
import logging
import snapatac2 as snap
import time
import multiprocessing as mp
import numpy as np
# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

start_time = time.time()
logging.info("Started: Feature QC")

# Inputs from Snakemake
try:
    adata_file = sys.argv[1]
    path_to_blacklist = sys.argv[2]
    filter_lower_quantile = float(sys.argv[3])
    filter_upper_quantile = float(sys.argv[4])
    n_features = int(sys.argv[5])
    output_file = sys.argv[6]
    output_features_txt = sys.argv[7]
except IndexError:
    logging.error("Not enough arguments provided.")
    sys.exit(1)

logging.info(f"Adata file: {adata_file}")
logging.info(f"Filter lower quantile: {filter_lower_quantile}")
logging.info(f"Filter upper quantile: {filter_upper_quantile}")
logging.info(f"N features: {n_features}")
logging.info(f"Output file: {output_file}")

# Load and copy to new output file
adata = snap.read(adata_file, backed='r', backend='hdf5')
adata_copy = adata.copy(filename=output_file, backend=None)
adata.close()

# 2. Initial Feature Selection (Required for Scrublet)
logging.info("Running initial feature selection for Scrublet...")
snap.pp.select_features(
    adata_copy,
    n_features=n_features,
    blacklist=path_to_blacklist,
    filter_lower_quantile=filter_lower_quantile,
    filter_upper_quantile=filter_upper_quantile,
    inplace=True,
    n_jobs=max(1, mp.cpu_count() - 2)
)

# 3. Doublet Detection & Removal
try:
    snap.pp.scrublet(adata_copy)  # Automatically uses .var['selected']
    snap.pp.filter_doublets(adata_copy)
    logging.info('Completed: Scrublet doublet removal.')
except Exception as e:
    logging.error(f'Scrublet failed with error: {e}')

# 4. Final Stratified Feature Selection (The "Pure" Vote)
# Recalculate variance now that the artificial doublets are gone
logging.info("Recalculating pure features on singlets only...")
snap.pp.select_features(
    adata_copy,
    n_features=n_features,
    blacklist=path_to_blacklist,
    max_iter=1,
    filter_lower_quantile=filter_lower_quantile,
    filter_upper_quantile=filter_upper_quantile,
    inplace=True,
    n_jobs=max(1, mp.cpu_count() - 2)
)


# 5. Extract and Save the "Vote"
# Force both the names and the mask into numpy arrays for safe boolean indexing
all_features = np.array(adata_copy.var_names)
mask = np.array(adata_copy.var['selected'])
sample_features = all_features[mask]

# sample_features = adata_copy.var_names[adata_copy.var['selected']].values

logging.info(f"Saving {len(sample_features)} features to {output_features_txt}...")
np.savetxt(output_features_txt, sample_features, fmt='%s')

adata_copy.close()

end_time = time.time()
duration = end_time - start_time
logging.info(f"Completed: Feature QC")
logging.info(f"Total running time: {duration:.2f} seconds")