import numpy as np
import pandas as pd
import sys
import os
import snapatac2 as snap
import multiprocessing as mp
import scanpy as sc
import logging
import time
import magic
# Ensure local modules can be imported
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    import python_utils  # Contains remove_genes
except ImportError:
    # Fallback or error if python_utils is missing
    logging.warning("Could not import python_utils. Ensure it is in the script directory.")

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

start_time = time.time()
logging.info("Started: Create GEM")

# --- Assign command-line arguments ---
try:
    h5ad_input = sys.argv[1]
    work_dir = sys.argv[2]
    genome_annot = sys.argv[3]
    min_cells = int(sys.argv[4])
    flavor = sys.argv[5]
    batch_key = sys.argv[6]
    n_top_genes = None if sys.argv[7].lower() == "none" else int(sys.argv[7])
    min_mean = float(sys.argv[8])
    max_mean = float(sys.argv[9])
    min_disp = float(sys.argv[10])
    n_jobs = int(sys.argv[11])
    organism = sys.argv[12]
    remove_mit = sys.argv[13].lower() == 'true'
    remove_ribo = sys.argv[14].lower() == 'true'
    remove_sex_genes = sys.argv[15].lower() == 'true'
    only_coding_genes = sys.argv[16].lower() == 'true'
    metadata_file_path = sys.argv[17]
    counts_output_h5ad = sys.argv[18]
except IndexError:
    logging.error("Insufficient arguments provided.")
    sys.exit(1)

logging.info(f"Dataset: {h5ad_input}")
logging.info(f"Work dir: {work_dir}")
logging.info(f"Genome annotation: {genome_annot}")
logging.info(f"Ensembl metadata: {metadata_file_path}")

# Change to work directory
if os.path.exists(work_dir):
    os.chdir(work_dir)
else:
    logging.warning(f"Work directory {work_dir} does not exist. Staying in current dir.")

# --- Create gene matrix ---
try:
    logging.info('Started: make_gene_matrix')

    # Load in memory (backed=None) as gene matrix creation requires full access
    data = snap.read_dataset(h5ad_input)

    # Generate Gene Matrix
    adata = snap.pp.make_gene_matrix(data, genome_annot)

    # Carry over UMAP from integration step
    if "X_umap_mnc" in data.obsm:
        adata.obsm["X_umap_mnc"] = data.obsm["X_umap_mnc"]
    if "X_umap_harmony" in data.obsm:
        adata.obsm["X_umap_harmony"] = data.obsm["X_umap_harmony"]
    if "X_umap" in data.obsm:
        adata.obsm["X_umap"] = data.obsm["X_umap"]
    data.close()

    # Modular Filtering
    if 'python_utils' in sys.modules:
        logging.info("Applying gene filtering using python_utils...")
        adata = python_utils.remove_genes(
            adata,
            organism=organism,
            remove_mit=remove_mit,
            remove_ribo=remove_ribo,
            remove_sex_genes=remove_sex_genes,
            only_coding_genes=only_coding_genes,
            metadata_path=metadata_file_path
        )
    else:
        logging.warning("Skipping detailed gene filtering (python_utils not loaded).")

    # Standard Scanpy filtering
    logging.info(f"Filtering features with min_cells={min_cells}...")
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.layers['counts'] = adata.X.copy()

    # --- Normalization ---
    logging.info("Normalizing and Logging...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers['normalized'] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()

    # --- Imputation ---
    logging.info('Started: MAGIC imputation')
    sc.external.pp.magic(adata, solver="approximate", n_jobs=n_jobs)
    logging.info('Finished: MAGIC')
    adata.layers['MAGIC'] = adata.X.copy()

    # Final cleanup: Remove spaces in var columns for H5AD compatibility
    adata.raw = None
    adata.var.columns = [c.replace(" ", "_").replace("/", "_") for c in adata.var.columns]

    logging.info(f"Writing AnnData to {counts_output_h5ad}...")
    adata.write_h5ad(counts_output_h5ad, compression="gzip")

except Exception as e:
    logging.error(f"Failed to create GEM: {e}")
    sys.exit(1)

end_time = time.time()
logging.info("Finished: Create GEM")
logging.info(f"Total running time: {end_time - start_time:.2f} seconds")