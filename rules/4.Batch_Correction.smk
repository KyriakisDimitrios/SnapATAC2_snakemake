# --- 1. Standard Branch ---
rule batch_std:
    input:
        flag = get_path("standard", "merge", "flag")
    params:
        dataset   = get_path("standard", "merge", "AnnDataSet"),
        blacklist = config['path_to_blacklist'],
        n_feat    = config['analysis_params']['batch']['n_features'],
        max_iter  = config['analysis_params']['batch']['max_iter'],
        n_iter    = config['analysis_params']['batch']['n_iter'],
        # Eigen Plot (Param because it is an output of the script but not a rule output used downstream)
        png_eigen = report(get_path("standard", "batch", "dir") + "eigen.png", category="Standard", subcategory="Batch")
    output:
        flag_out = get_path("standard", "batch", "flag"),
        png1     = report(get_path("standard", "batch", "dir") + "umap_sample.png", category="Standard", subcategory="Batch"),
        png2     = report(get_path("standard", "batch", "dir") + "umap_aft_sample.png", category="Standard", subcategory="Batch"),
        png3     = report(get_path("standard", "batch", "dir") + "umap_aft_leiden.png", category="Standard", subcategory="Batch")
    log:
        get_path("standard", "batch", "log")
    conda: '../envs/magic_env.yaml' # Assuming magic_env contains snapatac2/scanpy
    shell:
        """
        python scripts/4.Batch_Correction.py \
        {input.flag} \
        {params.dataset} \
        {params.blacklist} \
        {params.n_feat} \
        {params.max_iter} \
        {params.n_iter} \
        {params.png_eigen} \
        {output.png1} \
        {output.png2} \
        {output.png3} \
        {output.flag_out} \
        > {log} 2>&1
        """

# --- 2. Metadata Branch ---
rule batch_meta:
    input:
        flag = get_path("metadata", "merge", "flag")
    params:
        dataset   = get_path("metadata", "merge", "AnnDataSet"),
        blacklist = config['path_to_blacklist'],
        n_feat    = config['analysis_params']['batch']['n_features'],
        max_iter  = config['analysis_params']['batch']['max_iter'],
        n_iter    = config['analysis_params']['batch']['n_iter'],
        png_eigen = report(get_path("metadata", "batch", "dir") + "eigen.png", category="Metadata", subcategory="Batch")
    output:
        flag_out = get_path("metadata", "batch", "flag"),
        png1     = report(get_path("metadata", "batch", "dir") + "umap_sample.png", category="Metadata", subcategory="Batch"),
        png2     = report(get_path("metadata", "batch", "dir") + "umap_aft_sample.png", category="Metadata", subcategory="Batch"),
        png3     = report(get_path("metadata", "batch", "dir") + "umap_aft_leiden.png", category="Metadata", subcategory="Batch")
    log:
        get_path("metadata", "batch", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/4.Batch_Correction.py \
        {input.flag} \
        {params.dataset} \
        {params.blacklist} \
        {params.n_feat} \
        {params.max_iter} \
        {params.n_iter} \
        {params.png_eigen} \
        {output.png1} \
        {output.png2} \
        {output.png3} \
        {output.flag_out} \
        > {log} 2>&1
        """

# --- 3. DG Subset Branch ---
rule batch_DGsub:
    input:
        flag = get_path("DG_Subset", "merge", "flag")
    params:
        dataset   = get_path("DG_Subset", "merge", "AnnDataSet"),
        blacklist = config['path_to_blacklist'],
        n_feat    = config['analysis_params']['batch']['n_features'],
        max_iter  = config['analysis_params']['batch']['max_iter'],
        n_iter    = config['analysis_params']['batch']['n_iter'],
        png_eigen = report(get_path("DG_Subset", "batch", "dir") + "eigen.png", category="DG_Subset", subcategory="Batch")
    output:
        flag_out = get_path("DG_Subset", "batch", "flag"),
        png1     = report(get_path("DG_Subset", "batch", "dir") + "umap_sample.png", category="DG_Subset", subcategory="Batch"),
        png2     = report(get_path("DG_Subset", "batch", "dir") + "umap_aft_sample.png", category="DG_Subset", subcategory="Batch"),
        png3     = report(get_path("DG_Subset", "batch", "dir") + "umap_aft_leiden.png", category="DG_Subset", subcategory="Batch")
    log:
        get_path("DG_Subset", "batch", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/4.Batch_Correction.py \
        {input.flag} \
        {params.dataset} \
        {params.blacklist} \
        {params.n_feat} \
        {params.max_iter} \
        {params.n_iter} \
        {params.png_eigen} \
        {output.png1} \
        {output.png2} \
        {output.png3} \
        {output.flag_out} \
        > {log} 2>&1
        """