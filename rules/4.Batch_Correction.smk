# # # --- 1. Standard Branch ---
# rule batch_std:
#     input:
#         h5ad_input= get_path("standard","sample_qc_filter","h5ad_output")
#     params:
#         batch_var=config['analysis_params']['batch']['batch_var'],
#         blacklist=config['path_to_blacklist'],
#         n_feat=config['analysis_params']['batch']['n_features'],
#         max_iter=config['analysis_params']['batch']['max_iter'],
#         n_iter=config['analysis_params']['batch']['n_iter']
#     output:
#         flag_out=get_path("standard","batch","flag"),
#         h5ad_output=get_path("standard","batch","h5ad_output"),
#
#         # --- 1. Diagnostics ---
#         png_eigen=report(
#             get_path("standard","batch","dir") + "eigen.png",
#             category="Batch Correction",
#             subcategory="Diagnostics"
#         ),
#
#         # --- 2. Pre-Correction ---
#         png_pre=report(
#             get_path("standard","batch","dir") + "umap_pre_sample.png",
#             category="Batch Correction",
#             subcategory="Pre-Correction"
#         ),
#
#         # --- 3. MNC Method ---
#         png_mnc_s=report(
#             get_path("standard","batch","dir") + "umap_mnc_sample.png",
#             category="Batch Correction",
#             subcategory="MNC"
#         ),
#         png_mnc_l=report(
#             get_path("standard","batch","dir") + "umap_mnc_leiden.png",
#             category="Batch Correction",
#             subcategory="MNC"
#         ),
#
#         # --- 4. Harmony Method ---
#         png_harm_s=report(
#             get_path("standard","batch","dir") + "umap_harmony_sample.png",
#             category="Batch Correction",
#             subcategory="Harmony"
#         ),
#         png_harm_l=report(
#             get_path("standard","batch","dir") + "umap_harmony_leiden.png",
#             category="Batch Correction",
#             subcategory="Harmony"
#         )
#     log:
#         get_path("standard","batch","log")
#     conda: '../envs/snapatac2_env.yaml'
#     shell:
#         """
#         python scripts/4.Batch_Correction.py \
#         {input.h5ad_input} \
#         {params.batch_var} \
#         {params.blacklist} \
#         {params.n_feat} \
#         {params.max_iter} \
#         {params.n_iter} \
#         {output.png_eigen} \
#         {output.png_pre} \
#         {output.png_mnc_s} \
#         {output.png_mnc_l} \
#         {output.png_harm_s} \
#         {output.png_harm_l} \
#         {output.flag_out} \
#         {output.h5ad_output} \
#         > {log} 2>&1
#         """

#
# rule batch_iterative_std:
#     input:
#         h5ad_input= get_path("standard","sample_qc_filter","h5ad_output")
#     params:
#         batch_var=config['analysis_params']['batch']['batch_var'],
#         blacklist=config['path_to_blacklist'],
#         n_folds=config['analysis_params']['batch']['n_folds'],
#         cells_per_sample_fold=config['analysis_params']['batch']['cells_per_sample_fold'],
#         n_feat=config['analysis_params']['batch']['n_features'],
#         max_iter=config['analysis_params']['batch']['max_iter'],
#         n_iter=config['analysis_params']['batch']['n_iter']
#     output:
#         flag_out=get_path("standard","batch","flag"),
#         h5ad_output=get_path("standard","batch","h5ad_output"),
#
#         # --- 1. Diagnostics ---
#         png_eigen=report(
#             get_path("standard","batch","dir") + "eigen.png",
#             category="Batch Correction",
#             subcategory="Diagnostics"
#         ),
#
#         # --- 2. Pre-Correction ---
#         png_pre=report(
#             get_path("standard","batch","dir") + "umap_pre_sample.png",
#             category="Batch Correction",
#             subcategory="Pre-Correction"
#         ),
#
#         # --- 3. MNC Method ---
#         png_mnc_s=report(
#             get_path("standard","batch","dir") + "umap_mnc_sample.png",
#             category="Batch Correction",
#             subcategory="MNC"
#         ),
#         png_mnc_l=report(
#             get_path("standard","batch","dir") + "umap_mnc_leiden.png",
#             category="Batch Correction",
#             subcategory="MNC"
#         ),
#
#         # --- 4. Harmony Method ---
#         png_harm_s=report(
#             get_path("standard","batch","dir") + "umap_harmony_sample.png",
#             category="Batch Correction",
#             subcategory="Harmony"
#         ),
#         png_harm_l=report(
#             get_path("standard","batch","dir") + "umap_harmony_leiden.png",
#             category="Batch Correction",
#             subcategory="Harmony"
#         )
#     log:
#         get_path("standard","batch","log")
#     conda: '../envs/snapatac2_env.yaml'
#     shell:
#         """
#         python scripts/4.Iterative_Batch_Correction.py \
#         {input.h5ad_input} \
#         {params.batch_var} \
#         {params.blacklist} \
#         {params.n_folds} \
#         {params.cells_per_sample_fold} \
#         {params.n_feat} \
#         {params.max_iter} \
#         {params.n_iter} \
#         {output.png_eigen} \
#         {output.png_pre} \
#         {output.png_mnc_s} \
#         {output.png_mnc_l} \
#         {output.png_harm_s} \
#         {output.png_harm_l} \
#         {output.flag_out} \
#         {output.h5ad_output} \
#         > {log} 2>&1
#         """

rule batch_iterative_std:
    input:
        # Now correctly pulls the fully merged matrix from Rule 3
        h5ad_input = get_path("standard", "merge", "AnnDataSet")
    params:
        batch_var = config['analysis_params']['batch']['batch_var'],
        max_iter  = config['analysis_params']['batch']['max_iter'],
        n_iter    = config['analysis_params']['batch']['n_iter']
    output:
        flag_out    = get_path("standard", "batch", "flag"),
        h5ad_output = get_path("standard", "batch", "h5ad_output"),

        # --- 1. Diagnostics ---
        png_eigen = report(
            get_path("standard", "batch", "dir") + "eigen.png",
            category="Batch Correction",
            subcategory="Diagnostics"
        ),

        # --- 2. Pre-Correction ---
        png_pre = report(
            get_path("standard", "batch", "dir") + "umap_pre_sample.png",
            category="Batch Correction",
            subcategory="Pre-Correction"
        ),

        # --- 3. MNC Method ---
        png_mnc_s = report(
            get_path("standard", "batch", "dir") + "umap_mnc_sample.png",
            category="Batch Correction",
            subcategory="MNC"
        ),
        png_mnc_l = report(
            get_path("standard", "batch", "dir") + "umap_mnc_leiden.png",
            category="Batch Correction",
            subcategory="MNC"
        ),

        # --- 4. Harmony Method ---
        png_harm_s = report(
            get_path("standard", "batch", "dir") + "umap_harmony_sample.png",
            category="Batch Correction",
            subcategory="Harmony"
        ),
        png_harm_l = report(
            get_path("standard", "batch", "dir") + "umap_harmony_leiden.png",
            category="Batch Correction",
            subcategory="Harmony"
        )
    log:
        get_path("standard", "batch", "log")
    conda: '../envs/snapatac2_env.yaml'
    shell:
        """
        python scripts/4.Batch_Correction.py \
        {input.h5ad_input} \
        {params.batch_var} \
        {params.max_iter} \
        {params.n_iter} \
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