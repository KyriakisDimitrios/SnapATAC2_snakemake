import sys
import logging
import time
import warnings
import anndata as ad
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    stream=sys.stderr
)

# --- 0. Setup & Argument Parsing ---
scvi.settings.seed = 20120224
start_time = time.time()

try:
    args = sys.argv[1:]
    (reference_h5ad, query_h5ad, ref_celltype_col,
     min_cells, n_top_genes, n_layers, n_latent,
     max_epochs_scvi, max_epochs_scanvi, n_samples_per_label,
     png_scvi_elbo, png_scanvi_elbo, png_umap_integration,
     csv_output, flag_output) = args

    # Cast hyperparameters to integers
    min_cells = int(min_cells)
    n_top_genes = int(n_top_genes)
    n_layers = int(n_layers)
    n_latent = int(n_latent)
    max_epochs_scvi = int(max_epochs_scvi)
    max_epochs_scanvi = int(max_epochs_scanvi)
    n_samples_per_label = int(n_samples_per_label)

except Exception as e:
    logging.error(f"Argument Error: {e}")
    sys.exit(1)

logging.info(f"Using scvi-tools version: {scvi.__version__}")
logging.info(f"Started Integration and Annotation: {query_h5ad} vs {reference_h5ad}")

# --- 1. Load and Concatenate Data ---
reference = sc.read_h5ad(reference_h5ad)
query = sc.read_h5ad(query_h5ad)

if ref_celltype_col not in query.obs.columns:
    query.obs[ref_celltype_col] = pd.NA

data = ad.concat(
    [reference, query],
    join='inner',
    label='batch',
    keys=["reference", "query"],
    index_unique='_'
)

del reference
del query

# --- 2. Preprocessing ---
logging.info(f"Filtering (min_cells={min_cells}) and finding top {n_top_genes} variable genes...")
sc.pp.filter_genes(data, min_cells=min_cells)
sc.pp.highly_variable_genes(
    data,
    n_top_genes=n_top_genes,
    flavor="seurat_v3",
    batch_key="batch",
    subset=True
)

# --- 3. scVI Pre-training ---
logging.info("Setting up scVI AnnData...")
scvi.model.SCVI.setup_anndata(data, batch_key="batch")

logging.info(f"Training scVI VAE (layers={n_layers}, latent={n_latent}, epochs={max_epochs_scvi})...")
vae = scvi.model.SCVI(
    data,
    n_layers=n_layers,
    n_latent=n_latent,
    gene_likelihood="nb",
    dispersion="gene-batch",
)
vae.train(max_epochs=max_epochs_scvi, early_stopping=True)

# Save scVI ELBO Plot
fig, ax = plt.subplots(figsize=(8, 5))
vae.history['elbo_train'][1:].plot(ax=ax, label='Train ELBO')
if 'elbo_validation' in vae.history:
    vae.history['elbo_validation'].plot(ax=ax, label='Validation ELBO')
ax.set_title("scVI Training History")
ax.legend()
fig.savefig(png_scvi_elbo, bbox_inches='tight', dpi=300)
plt.close(fig)

# --- 4. scANVI Label Transfer ---
data.obs["celltype_scanvi"] = 'Unknown'
ref_idx = data.obs['batch'] == "reference"
data.obs.loc[ref_idx, "celltype_scanvi"] = data.obs.loc[ref_idx, ref_celltype_col]

logging.info(f"Training scANVI (epochs={max_epochs_scanvi}, samples/label={n_samples_per_label})...")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=data,
    labels_key="celltype_scanvi",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=max_epochs_scanvi, n_samples_per_label=n_samples_per_label)

# Save scANVI ELBO Plot
fig, ax = plt.subplots(figsize=(8, 5))
lvae.history['elbo_train'][1:].plot(ax=ax, label='Train ELBO')
ax.set_title("scANVI Training History")
ax.legend()
fig.savefig(png_scanvi_elbo, bbox_inches='tight', dpi=300)
plt.close(fig)

# --- 5. Prediction and Dimensionality Reduction ---
logging.info("Generating joint embedding...")
data.obs["predicted_cell_type"] = lvae.predict(data)
data.obsm["X_scANVI"] = lvae.get_latent_representation(data)

sc.pp.neighbors(data, use_rep="X_scANVI")
sc.tl.umap(data)

fig = sc.pl.umap(data, color=['predicted_cell_type', "batch"], wspace=0.45, show=False, return_fig=True)
fig.savefig(png_umap_integration, bbox_inches='tight', dpi=300)
plt.close(fig)

# --- 6. Extract Query Labels & Save ---
logging.info("Saving query annotations to CSV...")
query_obs = data.obs[data.obs['batch'] == 'query'].copy()

# Remove the '_query' suffix added by ad.concat to match original barcodes
query_obs.index = query_obs.index.str.replace('_query$', '', regex=True)

output_cols = ['predicted_cell_type']
query_obs[output_cols].to_csv(csv_output, index_label='barcode')

with open(flag_output, 'w') as f:
    f.write("done\n")

logging.info(f"Annotation Complete in {time.time() - start_time:.2f}s.")