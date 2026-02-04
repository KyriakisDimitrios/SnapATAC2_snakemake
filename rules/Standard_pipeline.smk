# ==============================================================================
# BRANCH A: STANDARD QC (Data Driven)
# Flow: Import -> Tile_QC -> Feature_QC -> Merge -> Batch -> Cluster -> GEM
# ==============================================================================

rule tile_qc_std:
    input:
        adatas = get_res_path(config['import_data']['h5ads_output_dir'] + "{sample}.h5ad"),
        annotation_gff3_file = config['genome_annot']
    output:
        tiled = get_res_path(config['standard']['tile_qc']['dir'] + "{sample}_Tiled.h5ad")
    params:
        min_tsse = config['analysis_params']['tile_qc']['min_tsse'],
        bin_size = config['analysis_params']['tile_qc']['bin_size']
    log: get_res_path(config['standard']['tile_qc']['log'] + "{sample}.log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/1.Tile_QC.py {input.adatas} {input.annotation_gff3_file} {params.min_tsse} {params.bin_size} {output.tiled} > {log} 2>&1"

rule feature_qc_std:
    input:
        adatas = get_res_path(config['standard']['tile_qc']['dir'] + "{sample}_Tiled.h5ad")
    params:
        low = config['analysis_params']['feature_qc']['filter_lower_quantile'],
        high = config['analysis_params']['feature_qc']['filter_upper_quantile'],
        n = config['analysis_params']['feature_qc']['n_features']
    output:
        features_out = get_res_path(config['standard']['feature_qc']['dir'] + "{sample}.h5ad")
    log: get_res_path(config['standard']['feature_qc']['log'] + "{sample}.log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/2.Feature_QC.py {input.adatas} {params.low} {params.high} {params.n} {output.features_out} > {log} 2>&1"


rule merge_std:
    input:
        expand(get_res_path(config['standard']['feature_qc']['dir'] + "{sample}.h5ad"), sample=SAMPLES),
        metadata_path = config['metadata_path']
    output:
        AnnDataSet = get_res_path(config['standard']['merge']['AnnDataSet']),
        flag = get_res_path(config['standard']['merge']['flag'])
    log: get_res_path(config['standard']['merge']['log'])
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/3.Merge_AnnData.py {input} {output.AnnDataSet} {output.flag} > {log} 2>&1"

rule batch_std:
    input: flag = get_res_path(config['standard']['merge']['flag'])
    params:
        dataset = get_res_path(config['standard']['merge']['AnnDataSet']),
        blacklist = config['path_to_blacklist'],
        n_feat = config['analysis_params']['batch']['n_features'],
        max_iter = config['analysis_params']['batch']['max_iter'],
        n_iter = config['analysis_params']['batch']['n_iter'],
        png = get_res_path(config['standard']['batch']['dir'] + "eigen.png")
    output:
        flag_out = get_res_path(config['standard']['batch']['flag']),
        png1 = get_res_path(config['standard']['batch']['dir'] + "umap_sample.png"),
        png2 = get_res_path(config['standard']['batch']['dir'] + "umap_aft_sample.png"),
        png3 = get_res_path(config['standard']['batch']['dir'] + "umap_aft_leiden.png")
    log: get_res_path(config['standard']['batch']['log'])
    conda: '../envs/scanpy_env.yaml'
    shell:
        "python scripts/4.Batch_Correction.py {input.flag} {params.dataset} {params.blacklist} {params.n_feat} {params.max_iter} {params.n_iter} {params.png} {output.png1} {output.png2} {output.png3} {output.flag_out} > {log} 2>&1"

rule cluster_std:
    input: flag = get_res_path(config['standard']['batch']['flag'])
    params:
        dataset = get_res_path(config['standard']['merge']['AnnDataSet']),
        res = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_res_path(config['standard']['clustering']['dir'])
    output:
        h5ad = get_res_path(config['standard']['clustering']['h5ad']),
        flag = get_res_path(config['standard']['clustering']['flag'])
    log: get_res_path(config['standard']['clustering']['log'])
    conda: '../envs/scanpy_env.yaml'
    shell:
        "python scripts/5.Clustering.py {input.flag} {params.dataset} {params.res} {params.out_dir} {output.h5ad} {output.flag} > {log} 2>&1"

rule gem_std:
    input: h5ad = get_res_path(config['standard']['clustering']['h5ad'])
    params:
        wd = config['projdir'],
        annot = config['genome_annot'],
        min_cells = config['analysis_params']['gem']['min_cells'],
        flavor = config['analysis_params']['gem']['flavor'],
        batch = config['analysis_params']['gem']['batch_key'],
        n_top = config['analysis_params']['gem']['n_top_genes'],
        min_m = config['analysis_params']['gem']['min_mean'],
        max_m = config['analysis_params']['gem']['max_mean'],
        min_d = config['analysis_params']['gem']['min_disp'],
        n_jobs = config['analysis_params']['gem']['n_jobs'],
        org = config['analysis_params']['gem']['organism'],
        rem_mit = config['analysis_params']['gem']['remove_mit'],
        rem_ribo = config['analysis_params']['gem']['remove_ribo'],
        rem_sex = config['analysis_params']['gem']['remove_sex_genes'],
        coding = config['analysis_params']['gem']['only_coding_genes'],
        meta = config['ensemble_database']
    output: gem = get_res_path(config['standard']['gem']['h5ad'])
    log: get_res_path(config['standard']['gem']['log'])
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/6.Create_GEM.py {input.h5ad} {params.wd} {params.annot} {params.min_cells} {params.flavor} {params.batch} {params.n_top} {params.min_m} {params.max_m} {params.min_d} {params.n_jobs} {params.org} {params.rem_mit} {params.rem_ribo} {params.rem_sex} {params.coding} {params.meta} {output.gem} > {log} 2>&1
        """

