import sys
import os
import warnings

import numpy as np
import pandas as pd

import scanpy as sc
import snapatac2 as snap

import multiprocessing as mp
import logging
import time
import concurrent.futures
from tqdm.auto import tqdm
# --- Logger Configuration ---
# Logs will appear instantly (unbuffered) in your Snakemake log file.
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    handlers=[logging.StreamHandler(sys.stderr)]
)


def load_data_to_memory(h5ad_input: str) -> snap.AnnData:
    """
    Materializes the full AnnDataSet into RAM.
    Run this cell ONCE to avoid repeating the 4-minute load time.
    """
    logging.info("1. Reading backed AnnDataSet pointers...")
    dat_backed = snap.read_dataset(h5ad_input)

    logging.info("2. Materializing full matrix into RAM (Paying the I/O tax once)...")
    dat_mem = dat_backed.to_adata()
    dat_backed.close()

    logging.info(f"Dataset fully in RAM. Shape: {dat_mem.shape}")
    return dat_mem


def split_to_memory_folds(
        dat_mem: snap.AnnData,
        batch_var: str = "sample",
        n_folds: int = 5,
        cells_per_sample: int = 2000
) -> list:
    """
    Rapidly slices the in-memory dataset into independent folds.
    Uses native logging to guarantee progress updates during heavy memory copying.
    """
    logging.info(f"3. Rapidly slicing {n_folds} folds...")
    batch_labels = dat_mem.obs[batch_var].to_numpy()
    unique_batches = np.unique(batch_labels)

    folds = []

    for fold_id in range(n_folds):
        start_fold = time.time()
        np.random.seed(fold_id * 42)
        keep_indices = []

        for b in unique_batches:
            idx = np.where(batch_labels == b)[0]
            if len(idx) > cells_per_sample:
                idx = np.random.choice(idx, cells_per_sample, replace=False)
            keep_indices.extend(idx)

        # The heavy C-level memory copy that freezes standard UI bars
        fold_adata = dat_mem[np.sort(keep_indices), :].copy()
        folds.append(fold_adata)

        # Bulletproof progress update
        elapsed = time.time() - start_fold
        logging.info(f"  ➔ Fold {fold_id + 1}/{n_folds} completed in {elapsed:.1f} seconds")

    logging.info("Successfully generated in-memory folds.")
    return folds



def _extract_fold_markers(
        fold_adata,  # snap.AnnData
        n_features,  # int
        sample_size,  # Union[int, None]
        n_jobs,  # int
        n_neighbors,  # int
        leiden_res,  # float
        n_top_markers  # int
):
    """
    Worker function: Processes a single fold dynamically.
    Using Union[int, None] for older Python compatibility.
    """
    logging.info(f"Processing fold with {fold_adata.shape[0]} cells...")

    snap.pp.select_features(fold_adata, n_features=n_features, inplace=True, n_jobs=n_jobs)
    snap.tl.spectral(fold_adata, features='selected', sample_size=sample_size, inplace=True)
    snap.pp.harmony(fold_adata, batch="sample", use_rep='X_spectral', key_added='X_spectral_harmony', n_jobs=n_jobs)
    snap.pp.knn(fold_adata, n_neighbors=n_neighbors, use_rep='X_spectral_harmony', inplace=True)
    snap.tl.leiden(fold_adata, key_added='leiden_pass1', resolution=leiden_res)

    # Extract Markers
    markers = snap.tl.marker_regions(fold_adata, groupby='leiden_pass1')

    # Functional map and union of markers
    extracted_features_map = map(lambda x: list(x)[:n_top_markers], markers.values())
    return reduce(lambda a, b: set(a).union(set(b)), extracted_features_map, set())


def mapped_robust_features_sequential(
        memory_folds,
        n_features=100000,
        sample_size=None,
        n_jobs=-1,
        n_neighbors=15,
        leiden_res=3.0,
        n_top_markers=2000
):
    """
    Maps the modular feature extraction pipeline sequentially across all folds.
    """
    logging.info(f"Sequentially mapping pipeline across {len(memory_folds)} folds...")

    worker = partial(
        _extract_fold_markers,
        n_features=n_features,
        sample_size=sample_size,
        n_jobs=n_jobs,
        n_neighbors=n_neighbors,
        leiden_res=leiden_res,
        n_top_markers=n_top_markers
    )

    fold_results = list(map(worker, memory_folds))
    final_features = reduce(lambda a, b: a.union(b), fold_results)

    logging.info(f"Mapping complete. Extracted {len(final_features)} robust features.")
    return np.array(list(final_features))


