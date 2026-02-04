rule features_qc:
    input:
        adatas = os.path.join(config['results_dir_path'], config['tile_qc']['tiled_output_dir'], "{sample}_Tiled.h5ad")
    params:
        filter_lower_quantile = config['features_qc']['filter_lower_quantile'],
        filter_upper_quantile = config['features_qc']['filter_upper_quantile'],
        n_features = config['features_qc']['n_features']
    output:
        features_out = os.path.join(config['results_dir_path'], config['features_qc']['h5ads_tiled_QC_snmk'], "{sample}.h5ad")
    log:
        log = os.path.join(config['results_dir_path'], config['features_qc']['log'], "2.features_qc_{sample}.log")
    shell:
        """
        source activate /sc/arion/projects/CommonMind/kyriad02/conda/envs/snapatac2Jan2025
        python workflow/scripts/2.Feature_QC.py \
            {input.adatas} \
            {params.filter_lower_quantile} \
            {params.filter_upper_quantile} \
            {params.n_features} \
            {output.features_out} > {log} 2>&1
        """