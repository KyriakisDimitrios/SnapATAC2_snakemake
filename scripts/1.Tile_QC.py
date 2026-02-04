import sys
import logging
import time
import snapatac2 as snap

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)

logging.info("Started: Tile QC")
start_time = time.time()

# Inputs from Snakemake
try:
    adata_file = sys.argv[1]
    annotation_gff3_file = sys.argv[2]
    min_tsse = int(sys.argv[3])
    bin_size = int(sys.argv[4])
    output_file = sys.argv[5]
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

# Compute TSS enrichment
snap.metrics.tsse(adata_copy, annotation_gff3_file)
logging.info('Completed: Compute TSS enrichment')

# Filter based on TSS
snap.pp.filter_cells(adata_copy, min_tsse=min_tsse)
logging.info('Completed: Filter based on TSS')

snap.pp.add_tile_matrix(adata_copy, bin_size=bin_size)
logging.info('Completed: Add Tile Matrix')

adata_copy.close()

end_time = time.time()
logging.info(f"Total time: {end_time - start_time:.2f}s")
logging.info("Completed: Tile QC")