# rule features_qc:
#     input:
#         adatas = os.path.join(config['results_dir_path'], config['tile_qc']['tiled_output_dir'], "{sample}_Tiled.h5ad")
#     params:
#         filter_lower_quantile = config['features_qc']['filter_lower_quantile'],
#         filter_upper_quantile = config['features_qc']['filter_upper_quantile'],
#         n_features = config['features_qc']['n_features']
#     output:
#         features_out = os.path.join(config['results_dir_path'], config['features_qc']['h5ads_tiled_QC_snmk'], "{sample}.h5ad")
#     log:
#         log = os.path.join(config['results_dir_path'], config['features_qc']['log'], "2.features_qc_{sample}.log")
#     shell:
#         """
#         source activate /sc/arion/projects/CommonMind/kyriad02/conda/envs/snapatac2Jan2025
#         python workflow/scripts/2.Feature_QC.py \
#             {input.adatas} \
#             {params.filter_lower_quantile} \
#             {params.filter_upper_quantile} \
#             {params.n_features} \
#             {output.features_out} > {log} 2>&1
#         """



# --- 2. Feature Selection (Standard Branch) ---
rule feature_qc_std:
    input:
        # Input is the output of the previous rule
        adatas = get_path("standard", "tile_qc", "h5ad_out")
    params:
        low  = config['analysis_params']['feature_qc']['filter_lower_quantile'],
        high = config['analysis_params']['feature_qc']['filter_upper_quantile'],
        n    = config['analysis_params']['feature_qc']['n_features']
    output:
        # Resolves to: results/standard/02.feature_qc/{sample}.h5ad
        # Note: We map this to the "feature_qc" structure definition
        features_out = os.path.join(get_path("standard", "feature_qc", "dir"), "{sample}.h5ad")
    log:
        get_path("standard", "feature_qc", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/2.Feature_QC.py \
        {input.adatas} \
        {params.low} \
        {params.high} \
        {params.n} \
        {output.features_out} \
        > {log} 2>&1
        """