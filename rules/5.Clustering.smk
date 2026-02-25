rule filter_batch_iterative_std:
    input:
        # Now correctly pulls the fully merged matrix from Rule 3
        h5ad_input=get_path("standard","batch","h5ad_output")
    params:
        batch_var=config['analysis_params']['filter_dr_cl']['batch_var'],
        n_features=config['analysis_params']['filter_dr_cl']['n_features'],
        max_iter=config['analysis_params']['filter_dr_cl']['max_iter'],
        n_iter=config['analysis_params']['filter_dr_cl']['n_iter'],
        exclude_clusters=config['analysis_params']['filter_dr_cl']['exclude_clusters']
    output:
        flag_out=get_path("standard","filter_dr_cl","flag"),
        h5ad_output=get_path("standard","filter_dr_cl","h5ad_output"),
        # csv_iterative= get_path("standard", "batch", "csv_iterative"),
        # --- 1. Diagnostics ---
        png_eigen=report(
            get_path("standard","filter_dr_cl","dir") + "eigen.png",
            category="Filtered Batch Correction",
            subcategory="Diagnostics"
        ),

        # --- 2. Pre-Correction ---
        png_pre=report(
            get_path("standard","filter_dr_cl","dir") + "umap_pre_sample.png",
            category="Filtered Batch Correction",
            subcategory="Pre-Correction"
        ),

        # --- 3. MNC Method ---
        png_mnc_s=report(
            get_path("standard","filter_dr_cl","dir") + "umap_mnc_sample.png",
            category="Filtered Batch Correction",
            subcategory="MNC"
        ),
        png_mnc_l=report(
            get_path("standard","filter_dr_cl","dir") + "umap_mnc_leiden.png",
            category="Filtered Batch Correction",
            subcategory="MNC"
        ),

        # --- 4. Harmony Method ---
        png_harm_s=report(
            get_path("standard","filter_dr_cl","dir") + "umap_harmony_sample.png",
            category="Filtered Batch Correction",
            subcategory="Harmony"
        ),
        png_harm_l=report(
            get_path("standard","filter_dr_cl","dir") + "umap_harmony_leiden.png",
            category="Filtered Batch Correction",
            subcategory="Harmony"
        )
    log:
        get_path("standard","filter_dr_cl","log")
    conda: '../envs/snapatac2_env.yaml'
    shell:
        """
        python scripts/5.Filter_DR_CL.py \
        {input.h5ad_input} \
        {params.batch_var} \
        {params.n_features} \
        {params.max_iter} \
        {params.n_iter} \
        {params.exclude_clusters} \
        {output.png_eigen} \
        {output.png_pre} \
        {output.png_mnc_s} \
        {output.png_mnc_l} \
        {output.png_harm_s} \
        {output.png_harm_l} \
        {output.flag_out} \
        {output.h5ad_output} \
        > {log} 2>&1
        """

#
# rule iterative_clustering_std:
#     input:
#         merged_h5ad=get_path("standard","merge","AnnDataSet")
#     output:
#         # Dynamically creates a separate folder for harmony vs mnc outputs
#         annotations_csv=os.path.join(config['results_dir_path'],"clustering","{algorithm}","iterative_annotations.csv")
#     params:
#         batch_key=config['analysis_params']['clustering']['batch_key'],
#         algorithm="{algorithm}"
#     log:
#         os.path.join(config['results_dir_path'],"Logs","iterative_clustering_{algorithm}.log")
#     conda: '../envs/snapatac2_env.yaml'
#     shell:
#         """
#         python scripts/5.Iterative_Clustering.py \
#         {input.merged_h5ad} \
#         {output.annotations_csv} \
#         {params.batch_key} \
#         {params.algorithm} \
#         > {log} 2>&1
#         """

#
# # --- 1. Standard Branch ---
# rule cluster_std:
#     input:
#         h5ad_input = get_path("standard","batch","h5ad_output"),
#     params:
#         res     = config['analysis_params']['clustering']['resolutions'],
#         out_dir = get_path("standard", "clustering", "dir")
#     output:
#         h5ad_output = get_path("standard", "clustering", "h5ad"),
#     log:
#         get_path("standard", "clustering", "log")
#     conda: '../envs/scanpy_env.yaml'
#     shell:
#         """
#         python scripts/5.Filter_DR_CL.py \
#         {input.h5ad_input} \
#         {params.res} \
#         {params.out_dir} \
#         {output.h5ad_output} > {log} 2>&1
#         """
#
# # --- 2. Metadata Branch ---
# rule cluster_meta:
#     input:
#         flag = get_path("metadata", "batch", "flag")
#     params:
#         dataset = get_path("metadata", "merge", "AnnDataSet"),
#         res     = config['analysis_params']['clustering']['resolutions'],
#         out_dir = get_path("metadata", "clustering", "dir")
#     output:
#         h5ad = get_path("metadata", "clustering", "h5ad"),
#         flag = get_path("metadata", "clustering", "flag")
#     log:
#         get_path("metadata", "clustering", "log")
#     conda: '../envs/scanpy_env.yaml'
#     shell:
#         """
#         python scripts/5.Filter_DR_CL.py \
#         {input.flag} \
#         {params.dataset} \
#         {params.res} \
#         {params.out_dir} \
#         {output.h5ad} \
#         {output.flag} \
#         > {log} 2>&1
#         """
#
# # --- 3. DG Subset Branch ---
# rule cluster_DGsub:
#     input:
#         flag = get_path("DG_Subset", "batch", "flag")
#     params:
#         dataset = get_path("DG_Subset", "merge", "AnnDataSet"),
#         res     = config['analysis_params']['clustering']['resolutions'],
#         out_dir = get_path("DG_Subset", "clustering", "dir")
#     output:
#         h5ad = get_path("DG_Subset", "clustering", "h5ad"),
#         flag = get_path("DG_Subset", "clustering", "flag")
#     log:
#         get_path("DG_Subset", "clustering", "log")
#     conda: '../envs/scanpy_env.yaml'
#     shell:
#         """
#         python scripts/5.Filter_DR_CL.py \
#         {input.flag} \
#         {params.dataset} \
#         {params.res} \
#         {params.out_dir} \
#         {output.h5ad} \
#         {output.flag} \
#         > {log} 2>&1
#         """