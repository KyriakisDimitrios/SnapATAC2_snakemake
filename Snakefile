# Description: Two-Branch Pipeline (Standard vs Metadata)
# Author: Dimitrios Kyriakis
# Date: January 30, 2026

import os
configfile: "config.yaml"

# --- Helper ---
def get_res_path(sub_dir):
    return os.path.join(config['results_dir_path'], sub_dir)

SAMPLES = config['samples']

# --- RULE 0: SHARED IMPORT ---
include: 'rules/0.import_data.smk'
include: 'rules/Standard_pipeline.smk'
include: 'rules/Constrain_data.smk'
include: 'rules/DG_subset.smk'
# --- DRIVER RULE ---
rule all:
    input:
        # get_res_path(config['standard']['gem']['h5ad']),
        # get_res_path(config['metadata']['gem']['h5ad']),
        get_res_path(config['DG_Subset']['peaks']['h5ad'])