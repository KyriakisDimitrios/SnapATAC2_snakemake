# Important note
# You read the object like this otherwise 0 vars (No peaks output)
# data = snap.read(h5ad_input, backed=None)
# snap.tl.macs3(data, groupby=groupby)
# peaks_df = snap.tl.merge_peaks(data.uns['macs3'], snap.genome.hg38)
# peak_mat = snap.pp.make_peak_matrix(data, use_rep=peaks_df['Peaks'],file=peaks_output_h5ad,backend='hdf5')

import sys
import os
import logging
import time
import snapatac2 as snap
import polars as pl
import warnings

# Multiprocessing guard remains essential
import multiprocessing


def main():
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
        logging.error("Insufficient arguments.")
        sys.exit(1)

    # Change to work directory
    if not os.path.exists(work_dir):
        os.makedirs(work_dir, exist_ok=True)
    os.chdir(work_dir)

    # --- 1. Load Data ---
    logging.info('Loading data...')
    data = snap.read(h5ad_input, backed=None)

    if groupby not in data.obs.columns:
        logging.error(f"CRITICAL: Column '{groupby}' not found in .obs!")
        sys.exit(1)

    # --- 2. Call Peaks (MACS3) ---
    logging.info(f'Calling peaks with MACS3 using {groupby}...')
    snap.tl.macs3(data, groupby=groupby)

    # --- 3. Merge Peaks ---
    logging.info('Merging peaks...')
    # Use hg38 or your specific genome
    peaks_df = snap.tl.merge_peaks(data.uns['macs3'], snap.genome.hg38)

    # --- NEW: Save Peaks to Disk ---
    # Construct a filename based on the output h5ad (e.g., results/peaks.tsv)
    peaks_tsv_path = peaks_output_h5ad.replace(".h5ad", "_merged_peaks.tsv")

    logging.info(f"Saving merged peaks to {peaks_tsv_path}...")
    try:
        if isinstance(peaks_df, pl.DataFrame):
            peaks_df.write_csv(peaks_tsv_path, separator='\t')
        else:
            # If it's a dict or list, convert to Polars first
            pl.DataFrame(peaks_df).write_csv(peaks_tsv_path, separator='\t')
    except Exception as e:
        logging.warning(f"Failed to save peaks CSV: {e}")

    # --- 4. Prepare List for Matrix ---
    if isinstance(peaks_df, pl.DataFrame):
        merged_peaks_list = peaks_df['Peaks'].to_list()
    elif isinstance(peaks_df, dict) and 'Peaks' in peaks_df:
        merged_peaks_list = peaks_df['Peaks']
    else:
        merged_peaks_list = peaks_df

    logging.info(f"Identified {len(merged_peaks_list)} merged peaks.")

    # --- 5. Create Peak Matrix ---
    logging.info('Creating cell-by-peak matrix...')
    peak_mat = snap.pp.make_peak_matrix(data, use_rep=peaks_df['Peaks'],file=peaks_output_h5ad,backend='hdf5')

    peak_mat.close()

    end_time = time.time()
    logging.info(f"Completed in {end_time - start_time:.2f}s")


# --- GUARD IS REQUIRED ---
if __name__ == "__main__":
    main()