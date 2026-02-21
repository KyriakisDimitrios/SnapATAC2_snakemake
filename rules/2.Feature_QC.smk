# # --- 2. Feature Selection (Standard Branch) ---
# rule feature_qc_std:
#     input:
#         # Input is the output of the previous rule
#         adatas = get_path("standard", "tile_qc", "h5ad_out"),
#
#     params:
#         low  = config['analysis_params']['feature_qc']['filter_lower_quantile'],
#         high = config['analysis_params']['feature_qc']['filter_upper_quantile'],
#         n    = config['analysis_params']['feature_qc']['n_features']
#     output:
#         # Resolves to: results/standard/02.feature_qc/{sample}.h5ad
#         # Note: We map this to the "feature_qc" structure definition
#         features_out = os.path.join(get_path("standard", "feature_qc", "dir"), "{sample}.h5ad")
#     log:
#         get_path("standard", "feature_qc", "log")
#     conda: '../envs/magic_env.yaml'
#     shell:
#         """
#         python scripts/2.Feature_QC.py \
#         {input.adatas} \
#         {params.low} \
#         {params.high} \
#         {params.n} \
#         {output.features_out} \
#         > {log} 2>&1
#         """

# --- 2. Feature Selection (Standard Branch) ---
rule feature_qc_std:
    input:
        # Input is the output of the previous rule
        adatas=get_path("standard","tile_qc","h5ad_out"),
        blacklist=config['path_to_blacklist']
    params:
        low=config['analysis_params']['feature_qc']['filter_lower_quantile'],
        high=config['analysis_params']['feature_qc']['filter_upper_quantile'],
        n=config['analysis_params']['feature_qc']['n_features']
    output:
        # 1. The filtered H5AD file
        features_out=os.path.join(get_path("standard","feature_qc","dir"),"{sample}.h5ad"),
        # 2. The new stratified feature "vote" text file
        features_txt=os.path.join(get_path("standard","feature_qc","dir"),"{sample}_features.txt")
    log:
        get_path("standard","feature_qc","log")
    conda: '../envs/snapatac2_env.yaml'
    shell:
        """
        python scripts/2.Feature_QC.py \
        {input.adatas} \
        {input.blacklist} \
        {params.low} \
        {params.high} \
        {params.n} \
        {output.features_out} \
        {output.features_txt} \
        > {log} 2>&1
        """
