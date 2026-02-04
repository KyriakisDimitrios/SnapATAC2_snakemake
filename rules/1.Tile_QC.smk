rule tile_qc:
    input:
        adatas = os.path.join(config['results_dir_path'], config['import_data']['h5ads_output_dir'], "{sample}.h5ad")
    params:
        annotation = config['genome_annot'],
        min_tsse = config['tile_qc']['min_tsse'],
        bin_size = config['tile_qc']['bin_size']
    output:
        tiled = os.path.join(config['results_dir_path'], config['tile_qc']['tiled_output_dir'], "{sample}_Tiled.h5ad")
    log:
        log = os.path.join(config['results_dir_path'], config['tile_qc']['log'], "1.tile_qc_{sample}.log")
    shell:
        """
        source activate /sc/arion/projects/CommonMind/kyriad02/conda/envs/snapatac2Jan2025
        python workflow/scripts/1.Tile_QC.py \
            {input.adatas} \
            {params.annotation} \
            {params.min_tsse} \
            {params.bin_size} \
            {output.tiled} > {log} 2>&1
        """