rule batch_QC:
    input:
        flag = os.path.join(config['results_dir_path'], config['convert_anndataset']['flag'])
    params:
        dataset = os.path.join(config['results_dir_path'], config['convert_anndataset']['AnnDataSet']),
        path_to_blacklist = config['path_to_blacklist'],
        n_features = config['batch_QC']["n_features"],
        max_iter = config['batch_QC']["max_iter"],
        n_iter = config['batch_QC']["n_iter"],
        png_eigenvalue = os.path.join(config['results_dir_path'], config['batch_QC']['png_eigenvalue'])
    output:
        png_sample = os.path.join(config['results_dir_path'], config['batch_QC']['png_sample']),
        png_aft_sample = os.path.join(config['results_dir_path'], config['batch_QC']['png_aft_sample']),
        png_aft_leiden = os.path.join(config['results_dir_path'], config['batch_QC']['png_aft_leiden']),
        flag_output = os.path.join(config['results_dir_path'], config['batch_QC']['flag_output'])
    log:
        log = os.path.join(config['results_dir_path'], config['batch_QC']['log'])
    shell:
        """
        source activate /sc/arion/projects/CommonMind/kyriad02/conda/envs/snapatac2Jan2025
        python workflow/scripts/4.Batch_Correction.py \
            {input.flag} \
            {params.dataset} \
            {params.path_to_blacklist} \
            {params.n_features} \
            {params.max_iter} \
            {params.n_iter} \
            {params.png_eigenvalue} \
            {output.png_sample} \
            {output.png_aft_sample} \
            {output.png_aft_leiden} \
            {output.flag_output} > {log} 2>&1
        """