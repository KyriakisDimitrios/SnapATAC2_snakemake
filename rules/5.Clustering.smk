rule iterative_clustering_std:
    input:
        merged_h5ad=get_path("standard","merge","AnnDataSet")
    output:
        # Dynamically creates a separate folder for harmony vs mnc outputs
        annotations_csv=os.path.join(config['results_dir_path'],"clustering","{algorithm}","iterative_annotations.csv")
    params:
        batch_key=config['analysis_params']['clustering']['batch_key'],
        algorithm="{algorithm}"
    log:
        os.path.join(config['results_dir_path'],"Logs","iterative_clustering_{algorithm}.log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/5.Iterative_Clustering.py \
        {input.merged_h5ad} \
        {output.annotations_csv} \
        {params.batch_key} \
        {params.algorithm} \
        > {log} 2>&1
        """


# --- 1. Standard Branch ---
rule cluster_std:
    input:
        h5ad_input = get_path("standard","batch","h5ad_output"),
    params:
        res     = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_path("standard", "clustering", "dir")
    output:
        h5ad_output = get_path("standard", "clustering", "h5ad"),
    log:
        get_path("standard", "clustering", "log")
    conda: '../envs/scanpy_env.yaml'
    shell:
        """
        python scripts/5.Clustering.py \
        {input.h5ad_input} \
        {params.res} \
        {params.out_dir} \
        {output.h5ad_output} > {log} 2>&1
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
    conda: '../envs/scanpy_env.yaml'
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
    conda: '../envs/scanpy_env.yaml'
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