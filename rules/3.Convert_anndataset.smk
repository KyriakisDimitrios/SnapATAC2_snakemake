rule convert_anndataset:
    input:
        processed_adatas = expand(
            os.path.join(config['results_dir_path'], config['features_qc']['h5ads_tiled_QC_snmk'], "{sample}.h5ad"),
            sample=SAMPLES
        ),
        metadata_path= config['metadata_path']
    output:
        AnnDataSet = os.path.join(config['results_dir_path'], config['convert_anndataset']['AnnDataSet']),
        flag = os.path.join(config['results_dir_path'], config['convert_anndataset']['flag'])
    log:
        log = os.path.join(config['results_dir_path'], config['convert_anndataset']['log'], "Merged_AnnData.log")
    shell:
        """
        source activate /sc/arion/projects/CommonMind/kyriad02/conda/envs/snapatac2Jan2025
        # Pass all inputs, then the two output paths
        python workflow/scripts/3.Merge_AnnData.py \
            {input.processed_adatas} \
            {input.metadata_path} \
            {output.AnnDataSet} \
            {output.flag} > {log} 2>&1
        """