import scanpy as sc
import numpy as np
from statsmodels.robust import mad  # Ensure this package is installed
import pandas as pd
import anndata
import warnings
import logging
import sys
import os
import scipy.io
import scanpy.external as sce

# --- Logger Configuration ---
# Logs will appear instantly (unbuffered) in your Snakemake log file.
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)


def load(path):
    """
    Load an AnnData object from an H5AD file and add sample metadata to the `obs` column.
    """
    logging.info(f"Loading AnnData from {path}...")
    try:
        # Read the H5AD file
        adata = sc.read_h5ad(path)
        # Extract sample metadata from the file path
        sample_name = path.split('/')[1].split('_')[0]
        adata.obs['Sample'] = sample_name
        logging.info(f"Loaded sample '{sample_name}' with shape: {adata.shape}")
        return adata
    except Exception as e:
        logging.error(f"Error loading {path}: {e}")
        raise


def load_sample_to_anndata(mtx_file_path):
    """
    Loads a single sample's gene score matrix and metadata into a Scanpy AnnData object.
    """
    logging.info(f"Processing {os.path.basename(mtx_file_path)}...")

    # Derive base name for metadata files
    base_name = os.path.splitext(os.path.basename(mtx_file_path))[0]
    cell_metadata_file = os.path.join(os.path.dirname(mtx_file_path), f"{base_name}_cell_metadata.tsv")
    gene_metadata_file = os.path.join(os.path.dirname(mtx_file_path), f"{base_name}_gene_metadata.tsv")

    # Check if all necessary files exist
    if not os.path.exists(cell_metadata_file):
        logging.warning(f"Cell metadata file not found for {base_name}: {cell_metadata_file}")
        return None
    if not os.path.exists(gene_metadata_file):
        logging.warning(f"Gene metadata file not found for {base_name}: {gene_metadata_file}")
        return None

    try:
        logging.info("Step 1: Reading the MTX data externally...")
        data_matrix = scipy.io.mmread(mtx_file_path)
        logging.info(f"Original matrix shape: {data_matrix.shape}")

        logging.info("Step 2: Transposing the matrix...")
        data_matrix_transposed = data_matrix.T
        logging.info(f"Transposed matrix shape: {data_matrix_transposed.shape}")

        cell_metadata = pd.read_csv(cell_metadata_file, sep="\t", index_col=0)
        logging.info(f"  Loaded cell metadata with shape: {cell_metadata.shape}")

        gene_metadata = pd.read_csv(gene_metadata_file, sep="\t", index_col=0)
        logging.info(f"  Loaded gene metadata with shape: {gene_metadata.shape}")

        if 'name' in gene_metadata.columns:
            gene_metadata = gene_metadata.set_index('name')

        adata = sc.AnnData(X=data_matrix_transposed,
                           obs=cell_metadata,
                           var=gene_metadata)

        adata.var_names_make_unique()
        adata.X = adata.X.tocsr()
        adata.layers['data'] = adata.X.copy()

        # Calculate standard QC metrics
        logging.info("Calculating initial QC metrics (Mito, Ribo, Hb)...")
        adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]|^Hb[^(P)]", regex=True)

        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
        )

        logging.info(f"  Successfully created AnnData object for {base_name}: {adata.shape}")
        return adata

    except Exception as e:
        logging.error(f"Error loading data for {base_name}: {e}")
        return None


def qc(adata, organism):
    """
    Add quality control (QC) metrics to an AnnData object.
    """
    logging.info(f"Running detailed QC for organism: {organism}")

    if organism == 'human':
        mt_var = "MT-"
        rpl_var = "RPL"
        rps_var = "RPS"
        hb_var = "HB-"
    else:
        mt_var = "Mt-"
        rpl_var = "Rpl-"
        rps_var = "Rps-"
        hb_var = "Hb-"

    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(mt_var)
    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith((rps_var, rpl_var))
    # Hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^" + hb_var + "[ab]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    logging.info("QC metrics calculated.")
    return adata


def remove_genes(adata, organism='human',
                 remove_mit=True,
                 remove_ribo=False,
                 remove_sex_genes=True,
                 only_coding_genes=True,
                 offline=True,
                 metadata_path=None):
    """
    Filter specific gene categories (mitochondrial, ribosomal, sex-linked, non-coding)
    from an AnnData object using a memory-efficient masking strategy.
    """
    # --- 1. Setup & Prefixes ---
    prefixes = {
        'human': {'mt': 'MT-', 'rpl': 'RPL', 'rps': 'RPS'},
        'mouse': {'mt': 'mt-', 'rpl': 'Rpl', 'rps': 'Rps'}
    }

    if organism not in prefixes:
        raise ValueError("Organism must be 'human' or 'mouse'")

    p = prefixes[organism]
    logging.info(f"Input Data: {adata.n_obs} cells x {adata.n_vars} genes")

    # --- 2. Metadata Handling ---
    if metadata_path is None:
        base_path = '/sc/arion/projects/CommonMind/kyriad02/Ensembl_Database/'
        metadata_path = f"{base_path}{organism}_ensembl_metadata.csv"

    meta_df = None

    try:
        if offline:
            # --- Offline Mode ---
            logging.info(f"Loading offline metadata from: {metadata_path}")
            meta_df = pd.read_csv(metadata_path)
        else:
            # --- Online Mode (BioMart) ---
            logging.info("Fetching metadata from BioMart...")
            from biomart import Dataset

            dataset_name = 'hsapiens_gene_ensembl' if organism == 'human' else 'mmusculus_gene_ensembl'
            dataset = Dataset(name=dataset_name, host='http://www.ensembl.org')

            meta_df = dataset.query(attributes=[
                'ensembl_gene_id', 'external_gene_name',
                'chromosome_name', 'gene_biotype'
            ])

            meta_df.rename(columns={
                'external_gene_name': 'Gene name',
                'ensembl_gene_id': 'Gene stable ID',
                'chromosome_name': 'Chromosome/scaffold name',
                'gene_biotype': 'Gene type'
            }, inplace=True)
            logging.info("BioMart fetch successful.")

        # --- Common Metadata Processing ---
        meta_df = meta_df.drop_duplicates(subset=['Gene name']).set_index('Gene name')
        adata.var = adata.var.join(meta_df, how='left')
        logging.info("Metadata successfully merged into adata.var")

    except Exception as e:
        logging.warning(f"Metadata error: {e}")
        logging.warning("Skipping metadata-dependent filters (Sex genes, Non-coding genes).")

    # --- 3. Create the "Master Mask" ---
    keep_mask = pd.Series(True, index=adata.var_names)

    # --- 4. Update the Mask (Logical AND) ---
    if remove_mit:
        is_mito = adata.var_names.str.startswith(p['mt'])
        keep_mask &= ~is_mito
        logging.info(f"Marked {sum(is_mito)} mitochondrial genes for removal.")

    if remove_ribo:
        is_ribo = adata.var_names.str.startswith((p['rpl'], p['rps']))
        keep_mask &= ~is_ribo
        logging.info(f"Marked {sum(is_ribo)} ribosomal genes for removal.")

    if remove_sex_genes:
        if 'Chromosome/scaffold name' in adata.var.columns:
            is_sex = adata.var['Chromosome/scaffold name'].isin(['X', 'Y'])
            keep_mask &= ~is_sex
            logging.info(f"Marked {sum(is_sex)} sex-linked genes for removal.")
        else:
            logging.warning("Skipping sex gene removal (Metadata column missing).")

    if only_coding_genes:
        if 'Gene type' in adata.var.columns:
            is_coding = adata.var['Gene type'] == 'protein_coding'
            keep_mask &= is_coding
            logging.info(f"Marked {sum(~is_coding)} non-coding/undefined genes for removal.")
        else:
            logging.warning("Skipping non-coding removal (Metadata column missing).")

    # --- 5. The Single Cut ---
    n_before = adata.n_vars
    adata_filtered = adata[:, keep_mask].copy()
    n_after = adata_filtered.n_vars

    logging.info(f"Filtering complete. Removed {n_before - n_after} genes.")
    logging.info(f"Final Output: {adata_filtered.n_obs} cells x {adata_filtered.n_vars} genes")

    return adata_filtered


def mad_outlier(adata, metric, nmads, upper_only=False):
    """
    Identify outliers in an AnnData object based on the Median Absolute Deviation (MAD).
    """
    M = adata.obs[metric]
    median_value = np.median(M)
    mad_value = mad(M)

    if not upper_only:
        return (M < median_value - nmads * mad_value) | (M > median_value + nmads * mad_value)
    return (M > median_value + nmads * mad_value)


def pp(adata, min_genes=200, nmad=5, mit_nmad=3, mit_thres=25):
    """
    Preprocess an AnnData object by removing outlier cells and filtering low-quality data.
    """
    logging.info("Starting Preprocessing (pp)...")
    logging.info(f"Filtering cells with < {min_genes} genes...")

    n_cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)

    logging.info("Identifying outliers using MAD...")
    bool_vector = (
            mad_outlier(adata, 'log1p_total_counts', nmad) +
            mad_outlier(adata, 'log1p_n_genes_by_counts', nmad) +
            mad_outlier(adata, 'pct_counts_in_top_20_genes', nmad) +
            mad_outlier(adata, 'pct_counts_mt', mit_nmad, upper_only=True)
    )

    adata.obs['Outlier'] = bool_vector
    logging.info(f"Marked {sum(bool_vector)} cells as outliers.")

    logging.info(f"Removing cells with > {mit_thres}% mitochondrial counts...")
    adata = adata[adata.obs.pct_counts_mt < mit_thres]

    n_cells_after = adata.n_obs
    logging.info(f"Preprocessing complete. Cells remaining: {n_cells_after} (Removed {n_cells_before - n_cells_after})")

    return adata


def filter_outliers(adata, var_name='Outlier'):
    logging.info(f"Total number of cells before outlier filtering: {adata.n_obs}")
    adata = adata[~adata.obs[var_name]]
    logging.info(f"Number of cells after filtering low quality cells: {adata.n_obs}")
    return adata


def dd(adata):
    """
    Identify and annotate potential doublets using Scrublet and DoubletDetection.
    """
    import doubletdetection

    logging.info("Starting Doublet Detection...")

    # Scrublet
    logging.info("Running Scrublet (expected_rate=0.1)...")
    sc.pp.scrublet(adata, expected_doublet_rate=0.1)

    # DoubletDetection
    logging.info("Running DoubletDetection BoostClassifier...")
    clf = doubletdetection.BoostClassifier(
        n_iters=10,
        clustering_algorithm="louvain",
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1
    )

    doublets = clf.fit(adata.X).predict(p_thresh=1e-3, voter_thresh=0.5)
    doublet_score = clf.doublet_score()

    adata.obs["clf_doublet"] = doublets
    adata.obs["clf_score"] = doublet_score

    logging.info("Doublet detection complete.")
    return adata


def rd(adata):
    """
    Remove identified potential doublets.
    """
    logging.info("Removing doublets...")
    n_before = adata.n_obs

    # Filter where either method identified a doublet
    adata = adata[~((adata.obs.predicted_doublet) | (adata.obs.clf_doublet == 1))]

    n_after = adata.n_obs
    logging.info(f"Removed {n_before - n_after} doublets. Cells remaining: {n_after}")
    return adata


def prep_int(adata):
    """
    Modify the `obs` index and save the raw counts in the AnnData object.
    """
    logging.info("Preparing data for integration...")

    # Make indices unique
    adata.obs.index = adata.obs.index + '_' + adata.obs.Sample

    # Save raw counts
    adata.layers['counts'] = adata.X.copy()

    logging.info("Unique indices created and raw counts saved to layers['counts'].")
    return adata





def get_mads(adata, column, mads=[4, 5], use_log10=False):
    """Calculates MAD-based thresholds on raw or log10 data."""
    try:
        data = pd.to_numeric(adata.obs[column], errors='coerce')
        # Filter out zeros/negatives for log transformation
        if use_log10:
            data = data[data > 0]
            data = np.log10(data)

        median_val = np.median(data)
        mad_val = np.median(np.abs(data - median_val))

        thresholds = {'median': median_val, 'mad': mad_val, 'is_log10': use_log10}

        for n in mads:
            thresholds[f'{n}_mad_upper'] = median_val + (n * mad_val)
            thresholds[f'{n}_mad_lower'] = median_val - (n * mad_val)

        return thresholds
    except Exception as e:
        logging.error(f"Error calculating MADs: {e}")
        return None


def save_qc_plot(adata, column, stats, use_log10, save_path, filename, plot_lower=False):
    """Generates and saves a histogram with optional lower threshold plotting."""
    try:
        data = pd.to_numeric(adata.obs[column], errors='coerce')

        if use_log10:
            data = data[data > 0]
            data = np.log10(data)
            label_suffix = " (log10)"
        else:
            label_suffix = ""

        plt.figure(figsize=(10, 6))
        plt.hist(data, bins=100, color='lightsteelblue', edgecolor='white', alpha=0.8)

        colors = ['red', 'darkred', 'orange', 'black']
        styles = ['--', '-.', ':', '-']

        # Get sorted MAD levels from the keys
        mad_levels = sorted(list(set([int(k.split('_')[0]) for k in stats.keys() if '_mad' in k])))

        for i, n in enumerate(mad_levels):
            # Upper threshold
            upper_val = stats[f'{n}_mad_upper']
            plt.axvline(upper_val, color=colors[i % len(colors)], linestyle=styles[i % len(styles)],
                        label=f'{n} MAD Upper ({upper_val:.2f})')

            # Lower threshold (only if requested)
            if plot_lower:
                lower_val = stats[f'{n}_mad_lower']
                # Use a dotted style for lower to distinguish it from upper
                plt.axvline(lower_val, color=colors[i % len(colors)], linestyle=':',
                            label=f'{n} MAD Lower ({lower_val:.2f})')

        plt.title(f"QC Distribution: {column}{label_suffix}", fontsize=14)
        plt.xlabel(f"{column}{label_suffix}", fontsize=12)
        plt.ylabel("Cell Count", fontsize=12)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  # Move legend outside if it gets too crowded
        plt.grid(axis='y', alpha=0.3)
        plt.tight_layout()

        os.makedirs(save_path, exist_ok=True)
        full_path = os.path.join(save_path, filename)
        plt.savefig(full_path, dpi=300)
        plt.close()
        logging.info(f"Plot successfully saved to: {full_path}")

    except Exception as e:
        logging.error(f"Error saving plot: {e}")

