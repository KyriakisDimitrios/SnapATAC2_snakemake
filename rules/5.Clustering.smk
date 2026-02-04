rule clustering:
    input:
        flag_input = os.path.join(config['results_dir_path'], config['batch_QC']['flag_output'])
    params:
        dataset = os.path.join(config['results_dir_path'], config['convert_anndataset']['AnnDataSet']),
        resolutions = config['clustering']['resolutions'],
        clustering_mnc = os.path.join(config['results_dir_path'], config['clustering']['clustering_mnc'])
    output:
        output_h5ad_path = os.path.join(config['results_dir_path'], config['clustering']['output_h5ad_path']),
        flag_output = os.path.join(config['results_dir_path'], config['clustering']['flag_output'])
    log:
        log = os.path.join(config['results_dir_path'], config['clustering']['log'])
    shell:
        """
        source activate /sc/arion/projects/CommonMind/kyriad02/conda/envs/snapatac2Jan2025
        python workflow/scripts/5.Clustering.py \
            {input.flag_input} \
            {params.dataset} \
            {params.resolutions} \
            {params.clustering_mnc} \
            {output.output_h5ad_path} \
            {output.flag_output} > {log} 2>&1
        """