import os
import sys
import logging
import time
import scanpy as sc
import snapatac2 as snap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from python_utils import get_mads
from python_utils import save_qc_plot

import warnings

# Suppress Deprecation warnings from background libraries
warnings.filterwarnings("ignore", category=DeprecationWarning)
# Suppress the "Divide by Zero" warning when doing log10
warnings.filterwarnings("ignore", message="divide by zero encountered in log10")

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

logging.info("Started: Tile QC")
start_time = time.time()

# --- Inputs from Snakemake ---
try:
    adata_file = sys.argv[1]
    annotation_gff3_file = sys.argv[2]
    path_to_blacklist = sys.argv[3]
    min_tsse = int(sys.argv[4])
    bin_size = int(sys.argv[5])
    output_file = sys.argv[6]
    output_qc_stats = sys.argv[7]
except IndexError:
    logging.error("Not enough arguments provided.")
    sys.exit(1)

logging.info(f"Adata file: {adata_file}")
logging.info(f"Annotation gff3 file: {annotation_gff3_file}")
logging.info(f"Min tsse: {min_tsse}")
logging.info(f"Bin size: {bin_size}")
logging.info(f"Output file: {output_file}")

adata = snap.read(adata_file, backed='r', backend='hdf5')
adata_copy = adata.copy(filename=output_file, backend=None)
adata.close()

logging.info(f"Starting Filtering. Original cells: {adata_copy.n_obs}")
snap.metrics.frip(adata_copy, {"blacklist_frac": path_to_blacklist}, inplace=True)

# Compute TSS enrichment
snap.metrics.tsse(adata_copy, annotation_gff3_file)
logging.info('Completed: Compute TSS enrichment')

# --- Extract Directories & Sample Name Early ---
base = os.path.basename(output_file)
sample_name = base.split('_')[0].replace('.h5ad', '')
directory_path = os.path.dirname(output_file)
qc_fig_dir = os.path.join(directory_path, "figures", "qc")

# Ensure the figure directory exists before we try to save plots to it
os.makedirs(qc_fig_dir, exist_ok=True)

# --- CAPTURE RAW METRICS (Before Filtering) ---
raw_ncells = adata_copy.n_obs
mean_tsse = adata_copy.obs['tsse'].mean()
median_tsse = adata_copy.obs['tsse'].median()
mean_nfrag = adata_copy.obs['n_fragment'].mean()
median_nfrag = adata_copy.obs['n_fragment'].median()

logging.info(f"Raw Metrics | Cells: {raw_ncells} | Median TSSE: {median_tsse:.2f} | Median Frags: {median_nfrag:.0f}")

# Save to CSV
qc_df = pd.DataFrame({
    'ncells': [raw_ncells],
    'mean_tsse': [mean_tsse],
    'median_tsse': [median_tsse],
    'mean_nfragments': [mean_nfrag],
    'median_nfragments': [median_nfrag]
}, index=[sample_name])

qc_df.to_csv(output_qc_stats, index_label='sample')
logging.info(f"Completed: Saved raw QC stats CSV to {output_qc_stats}")

# --- Plot 1: TSSe BEFORE Filtering ---
# before_png = os.path.join(qc_fig_dir, f"{sample_name}_tsse_scatter_before.png")
# snap.pl.tsse(adata_copy, interactive=False, out_file=before_png, height=500, width=500)
# logging.info(f"Completed: Saved pre-filtering TSSe plot to {before_png}")

# --- Plot 1: TSSe BEFORE Filtering ---
before_png = os.path.join(qc_fig_dir, f"{sample_name}_tsse_scatter_before.png")
# Bypass IPython by capturing the raw figure object and saving natively
fig_before = snap.pl.tsse(adata_copy, show=False)
fig_before.write_image(before_png, width=500, height=500)
logging.info(f"Completed: Saved pre-filtering TSSe plot to {before_png}")


# --- Proceed with Filtering ---
snap.pp.filter_cells(adata_copy, min_tsse=3)

# 1. Process TSSE (Raw scale)
tsse_stats = get_mads(adata=adata_copy,
                      column='tsse',
                      mads=[3, 4, 5],
                      use_log10=True)
save_qc_plot(adata=adata_copy,
             column='tsse',
             stats=tsse_stats,
             use_log10=True,
             save_path=qc_fig_dir,
             filename=sample_name + "_tsse_raw.png",
             plot_lower=True)

# 2. Process Fragments (Log10 scale)
frag_stats = get_mads(adata=adata_copy,
                      column='n_fragment',
                      mads=[3, 4, 5],
                      use_log10=True)
save_qc_plot(adata=adata_copy,
             column='n_fragment',
             stats=frag_stats,
             use_log10=True,
             save_path=qc_fig_dir,
             filename=sample_name + "_frag_log10.png",
             plot_lower=True)

# 3. Process Blacklist Fraction
blacklist_frac_stats = get_mads(adata=adata_copy,
                                column='blacklist_frac',
                                mads=[3, 4, 5],
                                use_log10=False)
save_qc_plot(adata=adata_copy,
             column='blacklist_frac',
             stats=blacklist_frac_stats,
             use_log10=False,
             save_path=qc_fig_dir,
             filename=sample_name + "_blacklist_frac.png",
             plot_lower=False)

# Log10 transformations for masking
tsse_data = pd.to_numeric(adata_copy.obs['tsse'], errors='coerce')
tsse_data = np.log10(tsse_data)
adata_copy.obs['log10(tsse)'] = tsse_data

n_fragment_data = pd.to_numeric(adata_copy.obs['n_fragment'], errors='coerce')
n_fragment_data = np.log10(n_fragment_data)
adata_copy.obs['log10(n_fragment)'] = n_fragment_data

# Apply MADs Mask
mask = (
        (adata_copy.obs['log10(tsse)'] > tsse_stats['4_mad_lower']) &
        (adata_copy.obs['log10(tsse)'] < tsse_stats['4_mad_upper']) &
        (adata_copy.obs['log10(n_fragment)'] < frag_stats['4_mad_upper']) &
        (adata_copy.obs['blacklist_frac'] < blacklist_frac_stats['5_mad_upper'])
)

adata_copy.subset(obs_indices=mask)
logging.info(f"MADs Filtering complete. Remaining cells: {adata_copy.n_obs}")

snap.pp.filter_cells(adata_copy, min_tsse=min_tsse)
logging.info('Completed: Final filter based on minimum TSSE')

# --- Plot 2: TSSe AFTER Filtering ---
# after_png = os.path.join(qc_fig_dir, f"{sample_name}_tsse_scatter_after.png")
# snap.pl.tsse(adata_copy, interactive=False, out_file=after_png, height=500, width=500)
# logging.info(f"Completed: Saved post-filtering TSSe plot to {after_png}")
# --- Plot 2: TSSe AFTER Filtering ---
after_png = os.path.join(qc_fig_dir, f"{sample_name}_tsse_scatter_after.png")
# Bypass IPython
fig_after = snap.pl.tsse(adata_copy, show=False)
fig_after.write_image(after_png, width=500, height=500)
logging.info(f"Completed: Saved post-filtering TSSe plot to {after_png}")

# --- Final Matrix Build ---
snap.pp.add_tile_matrix(adata_copy, bin_size=bin_size)
logging.info('Completed: Add Tile Matrix')

adata_copy.close()

end_time = time.time()
logging.info(f"Total time: {end_time - start_time:.2f}s")
logging.info("Completed: Tile QC")

# import os
# import sys
# import logging
# import time
# import scanpy as sc
# import snapatac2 as snap
#
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
#
# from python_utils import get_mads
# from python_utils import save_qc_plot
#
# import warnings
# # Suppress Deprecation warnings from background libraries
# warnings.filterwarnings("ignore", category=DeprecationWarning)
# # Suppress the "Divide by Zero" warning when doing log10
# warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
#
#
#
# # --- Logger Configuration ---
# logging.basicConfig(
#     format='%(asctime)s | %(levelname)-7s | %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     level=logging.INFO,
#     handlers=[logging.StreamHandler(sys.stderr)]
# )
#
# logging.info("Started: Tile QC")
# start_time = time.time()
#
# # Inputs from Snakemake
# try:
#     adata_file = sys.argv[1]
#     annotation_gff3_file = sys.argv[2]
#     path_to_blacklist = sys.argv[3]
#     min_tsse = int(sys.argv[4])
#     bin_size = int(sys.argv[5])
#     output_file = sys.argv[6]
#     output_qc_stats = sys.argv[7]
# except IndexError:
#     logging.error("Not enough arguments provided.")
#     sys.exit(1)
#
# logging.info(f"Adata file: {adata_file}")
# logging.info(f"Annotation gff3 file: {annotation_gff3_file}")
# logging.info(f"Min tsse: {min_tsse}")
# logging.info(f"Bin size: {bin_size}")
# logging.info(f"Output file: {output_file}")
#
# adata = snap.read(adata_file, backed='r', backend='hdf5')
# adata_copy = adata.copy(filename=output_file, backend=None)
# adata.close()
#
#
# logging.info(f"Starting Filtering. Original cells: {adata_copy.n_obs}")
# snap.metrics.frip(adata_copy, {"blacklist_frac": path_to_blacklist},inplace=True)
#
# # Compute TSS enrichment
# snap.metrics.tsse(adata_copy, annotation_gff3_file)
# logging.info('Completed: Compute TSS enrichment')
#
# # 2. CAPTURE RAW METRIC (Before Filtering)
# raw_mean_tsse = adata_copy.obs['tsse'].mean()
# raw_cell_count = adata_copy.n_obs
#
# logging.info(f"Raw Mean TSSE: {raw_mean_tsse:.4f}")
#
# # Save simple text file: "SCORE,CELL_COUNT"
# with open(output_qc_stats, 'w') as f:
#     f.write(f"{raw_mean_tsse},{raw_cell_count}")
#
# snap.pp.filter_cells(adata_copy, min_tsse=3)
#
# # 1. Get the filename: 'Rush-01_Tiled.h5ad'
# base = os.path.basename(output_file)
# # 2. Split by the underscore and take the first part: 'Rush-01'
# sample_name = base.split('_')[0]
#
# # Extract the directory path
# directory_path = os.path.dirname(output_file)
#
# # 1. Process TSSE (Raw scale)
# tsse_stats = get_mads(adata=adata_copy,
#          column='tsse',
#          mads=[3,4,5],
#          use_log10=True)
# save_qc_plot(adata=adata_copy,
#              column='tsse',
#              stats=tsse_stats,
#              use_log10=True,
#              save_path=directory_path+"/figures/qc",
#              filename=sample_name+"_tsse_raw.png",
#              plot_lower=True)
#
# # 2. Process Fragments (Log10 scale)
# frag_stats = get_mads(adata=adata_copy,
#          column='n_fragment',
#          mads=[3,4, 5],
#          use_log10=True)
# save_qc_plot(adata=adata_copy,
#              column='n_fragment',
#              stats=frag_stats,
#              use_log10=True,
#              save_path=directory_path+"/figures/qc",
#              filename=sample_name+"_frag_log10.png",
#              plot_lower=True)
#
# blacklist_frac_stats = get_mads(adata=adata_copy,
#          column='blacklist_frac',
#          mads=[3,4,5],
#          use_log10=False)
# save_qc_plot(adata=adata_copy,
#              column='blacklist_frac',
#              stats=blacklist_frac_stats,
#              use_log10=False,
#              save_path=directory_path+"/figures/qc",
#              filename=sample_name+"_blacklist_frac.png",
#              plot_lower=False)
#
# tsse_data = pd.to_numeric(adata_copy.obs['tsse'], errors='coerce')
# tsse_data = np.log10(tsse_data)
# adata_copy.obs['log10(tsse)'] = tsse_data
# n_fragment_data = pd.to_numeric(adata_copy.obs['n_fragment'], errors='coerce')
# n_fragment_data = np.log10(n_fragment_data)
# adata_copy.obs['log10(n_fragment)'] = n_fragment_data
#
#
# # Assuming tsse_stats and frag_stats were calculated with use_log10=True
# import logging
#
# # Assuming tsse_stats and frag_stats were calculated with use_log10=True
# mask = (
#     (adata_copy.obs['log10(tsse)'] > tsse_stats['4_mad_lower']) &
#     (adata_copy.obs['log10(tsse)'] < tsse_stats['4_mad_upper']) &
#     (adata_copy.obs['log10(n_fragment)'] < frag_stats['4_mad_upper'])&
#     (adata_copy.obs['blacklist_frac'] < blacklist_frac_stats['5_mad_upper'])
# )
#
# adata_copy.subset(obs_indices=mask)
# logging.info(f"Filtering complete. Remaining cells: {adata_copy.n_obs}")
#
# snap.pp.filter_cells(adata_copy, min_tsse=min_tsse)
# logging.info('Completed: Filter based on TSS')
#
# snap.pp.add_tile_matrix(adata_copy, bin_size=bin_size)
# logging.info('Completed: Add Tile Matrix')
#
# adata_copy.close()
#
# end_time = time.time()
# logging.info(f"Total time: {end_time - start_time:.2f}s")
# logging.info("Completed: Tile QC")