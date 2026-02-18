# --- 1. Tile & QC (Standard Branch) ---
rule tile_qc_std:
    input:
        # FIX: Must match rule import_data output exactly
        adatas = os.path.join(config['results_dir_path'], config['import_data']['h5ads_output_dir'], "{sample}.h5ad"),
        annotation_gff3_file = config['genome_annot'],
        path_to_blacklist= config['path_to_blacklist']
    params:
        min_tsse=config['analysis_params']['tile_qc']['min_tsse'],
        bin_size=config['analysis_params']['tile_qc']['bin_size']
    output:
        tiled = get_path("standard", "tile_qc", "h5ad_out"),
        raw_stats= get_path("standard", "tile_qc", "raw_qc")
    log:
        get_path("standard", "tile_qc", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/1.Tile_QC.py \
        {input.adatas} \
        {input.annotation_gff3_file} \
        {input.path_to_blacklist} \
        {params.min_tsse} \
        {params.bin_size} \
        {output.tiled} \
        {output.raw_stats} \
        > {log} 2>&1
        """
#
# # --- 1. Tile & QC (Standard Branch) ---
# rule tile_qc_std:
#     input:
#         # Input comes from the shared "import_data" step (00.raw_fragments)
#         # Note: We construct this path manually since 'import_data' isn't in the 'structure' block yet,
#         # or we can just use the config path directly since it's shared.
#         adatas = os.path.join(config['projdir'], config['import_data']['h5ads_output'].format(sample="{sample}")),
#         annotation_gff3_file = config['genome_annot']
#     output:
#         # Resolves to: results/standard/01.tiled/{sample}_Tiled.h5ad
#         tiled = get_path("standard", "tile_qc", "h5ad_out")
#     params:
#         min_tsse = config['analysis_params']['tile_qc']['min_tsse'],
#         bin_size = config['analysis_params']['tile_qc']['bin_size']
#     log:
#         get_path("standard", "tile_qc", "log")
#     conda: '../envs/magic_env.yaml'
#     shell:
#         """
#         python scripts/1.Tile_QC.py \
#         {input.adatas} \
#         {input.annotation_gff3_file} \
#         {params.min_tsse} \
#         {params.bin_size} \
#         {output.tiled} \
#         > {log} 2>&1
#         """
#
