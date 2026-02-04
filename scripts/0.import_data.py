import os
import sys
import logging
import snapatac2 as snap
import pandas as pd
import scanpy as sc
import time


def setup_logging():
    """Configures logging to write to stderr: timestamp | level | message."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        stream=sys.stderr
    )


def run_import():
    logger = logging.getLogger(__name__)
    print(sys.argv)
    # Inputs from Snakemake positional arguments
    try:
        frags = sys.argv[1]
        metadata_path = sys.argv[2]
        sample_name = sys.argv[3]
        min_frags = int(sys.argv[4])
        n_jobs = int(sys.argv[5])
        output_file = sys.argv[6]
    except IndexError:
        logger.error("Missing arguments. Expected: frags, metadata_path, sample_name, min_frags, n_jobs, output_file")
        sys.exit(1)

    logger.info("Started: Import Fragments")
    logger.info(f"Fragment file: {frags}")
    logger.info(f"Min num fragments: {min_frags}")
    logger.info(f"Output file: {output_file}")

    try:
        logging.info("Import fragments...")
        adata = snap.pp.import_fragments(
            frags,
            file=output_file,
            chrom_sizes=snap.genome.hg38,
            min_num_fragments=min_frags,
            sorted_by_barcode=False,
            n_jobs=n_jobs
        )
        logging.info("Merging metadata columns...")
        adata.obs_names = [f"{bc}_{sample_name}" for bc in adata.obs_names]
        meta_df = pd.read_csv(metadata_path, index_col='Unnamed: 0')
        meta_df['CellID'] = meta_df.index
        meta_subset = meta_df.reindex(adata.obs_names)
        for col in meta_subset.columns:
            series = meta_subset[col]
            if series.dtype == 'object' or series.dtype.name == 'category':
                values = series.astype(str).values
            else:
                values = series.values
            adata.obs[col] = values

        # The join happens on the Index automatically
        #adata.obs = adata.obs.join(meta_df, how='left')

        adata.close()



        # adata = sc.read_h5ad(output_file)
        # logging.info(f"Renaming observations with suffix _{sample_name}...")
        # adata.obs_names = [f"{bc}_{sample_name}" for bc in adata.obs_names]
        #
        # # --- 2. Load Metadata ---
        # logging.info(f"Loading metadata from {metadata_path}...")
        # meta_df = pd.read_csv(metadata_path, index_col='Unnamed: 0')
        # meta_df['CellID'] = meta_df.index
        # if ':' in adata.obs_names[0]:
        #     logging.info("Detected merge prefix in adata index. Cleaning...")
        #     # Split by ':' and take the second part (the real barcode_sample)
        #     clean_index = adata.obs_names.to_series().str.split(':').str[1]
        #     adata.obs_names = clean_index
        #
        # # 4. Merge
        # logging.info("Merging metadata columns...")
        # # The join happens on the Index automatically
        # adata.obs = adata.obs.join(meta_df, how='left')
        # adata.write(output_file)

        logger.info("Completed: Import Fragments")
    except Exception as e:
        logger.error(f"Failed to import fragments: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    setup_logging()
    run_import()