# --- 1. Standard Branch ---
rule merge_std:
    input:
        adatas = expand(
            get_path("standard", "feature_qc", "dir") + "{sample}.h5ad",
            sample=config['samples']
        ),
        meta_csv = config['metadata_path']
    output:
        AnnDataSet = get_path("standard", "merge", "AnnDataSet"),
        flag       = get_path("standard", "merge", "flag")
    log:
        get_path("standard", "merge", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        # ORDER MUST MATCH PYTHON UNPACKING:
        # *processed_adatas (All inputs), metadata, output, flag
        "python scripts/3.Merge_AnnData.py {input.adatas} {input.meta_csv} {output.AnnDataSet} {output.flag} > {log} 2>&1"


# --- 2. Metadata Branch ---
rule merge_meta:
    input:
        adatas = expand(
            os.path.join(config['results_dir'], config['branches']['metadata']['constrain']['dir'], "{sample}.h5ad"),
            sample=config['samples']
        ),
        meta_csv = config['metadata_path']
    output:
        AnnDataSet = get_path("metadata", "merge", "AnnDataSet"),
        flag       = get_path("metadata", "merge", "flag")
    log:
        get_path("metadata", "merge", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/3.Merge_AnnData.py {input.adatas} {input.meta_csv} {output.AnnDataSet} {output.flag} > {log} 2>&1"


# --- 3. DG_Subset Branch ---
rule merge_DGsub:
    input:
        adatas = expand(
            get_path("DG_Subset", "subset", "h5ad"),
            sample=config['samples']
        ),
        meta_csv = config['metadata_path']
    output:
        AnnDataSet = get_path("DG_Subset", "merge", "AnnDataSet"),
        flag       = get_path("DG_Subset", "merge", "flag")
    log:
        get_path("DG_Subset", "merge", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/3.Merge_AnnData.py {input.adatas} {input.meta_csv} {output.AnnDataSet} {output.flag} > {log} 2>&1"