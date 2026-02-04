import sys
import logging
import snapatac2 as snap
import time

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
    filter_lower_quantile = float(sys.argv[2])
    filter_upper_quantile = float(sys.argv[3])
    n_features = int(sys.argv[4])
    output_file = sys.argv[5]
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

# Select Features
snap.pp.select_features(
    adata_copy,
    filter_lower_quantile=filter_lower_quantile,
    filter_upper_quantile=filter_upper_quantile,
    n_features=n_features
)
logging.info('Completed: Select features')

# Run Scrublet for doublet detection
try:
    snap.pp.scrublet(adata_copy)
    logging.info('Completed: Run Scrublet for doublet detection.')
except Exception as e:
    logging.error(f'Scrublet failed with error: {e}')

# Remove detected doublets
snap.pp.filter_doublets(adata_copy)
logging.info('Completed: Remove detected doublets from the dataset.')

adata_copy.close()

end_time = time.time()
duration = end_time - start_time
logging.info(f"Completed: Feature QC")
logging.info(f"Total running time: {duration:.2f} seconds")