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
        sample_name = sys.argv[2]
        min_frags = int(sys.argv[3])
        n_jobs = int(sys.argv[4])
        output_file = sys.argv[5]
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

        adata.close()
        logger.info("Completed: Import Fragments")
    except Exception as e:
        logger.error(f"Failed to import fragments: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    setup_logging()
    run_import()