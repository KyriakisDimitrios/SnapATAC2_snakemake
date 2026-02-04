import sys
import os
import logging
import time
import snapatac2 as snap
import polars as pl

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

start_time = time.time()
logging.info(f"Started: Create Peaks (SnapATAC2 v{snap.__version__})")

# --- Assign command-line arguments ---
try:
    h5ad_input = sys.argv[1]
    work_dir = sys.argv[2]
    groupby = sys.argv[3]
    peaks_output_h5ad = sys.argv[4]
except IndexError:
    logging.error("Insufficient arguments. Usage: python create_peaks.py <input> <work_dir> <groupby> <output>")
    sys.exit(1)

logging.info(f"Dataset: {h5ad_input}")
logging.info(f"Grouping by: {groupby}")
logging.info(f"Output: {peaks_output_h5ad}")

# Change to work directory (useful for MACS3 temp files)
if os.path.exists(work_dir):
    os.chdir(work_dir)

# --- 1. Load Data ---
# Loading in memory (backed=None) is required for MACS3 in some versions,
# but for very large datasets, ensure you have enough RAM.
logging.info('Loading data...')
data = snap.read(h5ad_input, backed=None)

# --- 2. Call Peaks (MACS3) ---
logging.info('Calling peaks with MACS3...')
# This stores results in data.uns['macs3']
snap.tl.macs3(data, groupby=groupby)

# --- 3. Merge Peaks ---
logging.info('Merging peaks...')
# snap.genome.hg38 provides the chromosome sizes for human.
# Ensure this matches your organism.
peaks_df = snap.tl.merge_peaks(data.uns['macs3'], snap.genome.hg38)

# The output is a Polars DataFrame or Dictionary. We need the actual peak list.
# Usually 'Peaks' column holds the unique merged regions.
if isinstance(peaks_df, pl.DataFrame):
    # Convert Polars column to list
    merged_peaks_list = peaks_df['Peaks'].to_list()
elif isinstance(peaks_df, dict) and 'Peaks' in peaks_df:
    merged_peaks_list = peaks_df['Peaks']
else:
    # Fallback if it returns a list directly
    merged_peaks_list = peaks_df

logging.info(f"Identified {len(merged_peaks_list)} merged peaks.")

# --- 4. Create Peak Matrix ---
logging.info('Creating cell-by-peak matrix...')

peak_mat = snap.pp.make_peak_matrix(
    data,
    peaks=merged_peaks_list,
    file=peaks_output_h5ad,  # Directly save to the output path
    backend='hdf5'           # Use HDF5 backend
)

# --- 5. Clean up and Close ---
# peak_mat is now backed on disk at 'peaks_output_h5ad'
peak_mat.close()
data.close()

end_time = time.time()
logging.info(f"Completed in {end_time - start_time:.2f}s")