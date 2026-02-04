# Description: Two-Branch Pipeline (Standard vs Metadata)
# Author: Dimitrios Kyriakis
# Date: January 30, 2026

import os
configfile: "BICCN_config.yaml"


import os

# 1. Load the Public Config (Structure & Logic)
configfile: "config.yaml"

# 2. Load the Private Config (Sensitive Paths) if it exists
if os.path.exists("config_private.yaml"):
    configfile: "config_private.yaml"
else:
    # Optional: Warn if running locally without private data
    print("Warning: 'config_private.yaml' not found. Using defaults/placeholders.")



# --- Helper ---
def get_res_path(sub_dir):
    return os.path.join(config['results_dir_path'], sub_dir)

def get_path(branch_name, step, key):
    root = config['branches'][branch_name]['root']
    rel_path = config['structure'][step][key]
    return os.path.join(config['results_dir'], root, rel_path)

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