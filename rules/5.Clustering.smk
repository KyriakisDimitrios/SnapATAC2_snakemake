# --- 1. Standard Branch ---
rule cluster_std:
    input:
        flag = get_path("standard", "batch", "flag")
    params:
        # Note: Clustering reads the merged file (which was modified in-place by Batch Correction)
        dataset = get_path("standard", "merge", "AnnDataSet"),
        res     = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_path("standard", "clustering", "dir")
    output:
        h5ad = get_path("standard", "clustering", "h5ad"),
        flag = get_path("standard", "clustering", "flag")
    log:
        get_path("standard", "clustering", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/5.Clustering.py \
        {input.flag} \
        {params.dataset} \
        {params.res} \
        {params.out_dir} \
        {output.h5ad} \
        {output.flag} \
        > {log} 2>&1
        """

# --- 2. Metadata Branch ---
rule cluster_meta:
    input:
        flag = get_path("metadata", "batch", "flag")
    params:
        dataset = get_path("metadata", "merge", "AnnDataSet"),
        res     = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_path("metadata", "clustering", "dir")
    output:
        h5ad = get_path("metadata", "clustering", "h5ad"),
        flag = get_path("metadata", "clustering", "flag")
    log:
        get_path("metadata", "clustering", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/5.Clustering.py \
        {input.flag} \
        {params.dataset} \
        {params.res} \
        {params.out_dir} \
        {output.h5ad} \
        {output.flag} \
        > {log} 2>&1
        """

# --- 3. DG Subset Branch ---
rule cluster_DGsub:
    input:
        flag = get_path("DG_Subset", "batch", "flag")
    params:
        dataset = get_path("DG_Subset", "merge", "AnnDataSet"),
        res     = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_path("DG_Subset", "clustering", "dir")
    output:
        h5ad = get_path("DG_Subset", "clustering", "h5ad"),
        flag = get_path("DG_Subset", "clustering", "flag")
    log:
        get_path("DG_Subset", "clustering", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/5.Clustering.py \
        {input.flag} \
        {params.dataset} \
        {params.res} \
        {params.out_dir} \
        {output.h5ad} \
        {output.flag} \
        > {log} 2>&1
        """