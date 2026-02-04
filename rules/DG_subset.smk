# ==============================================================================
# BRANCH B: METADATA CONSTRAINED (Whitelist)
# Flow: Import -> Constrain -> Merge -> Batch -> Cluster -> GEM
# ==============================================================================
rule subset_DGsub:
    input:
        adata = get_res_path(config['import_data']['h5ads_output_dir'] + "{sample}.h5ad"),
    params:
        annotation = config['genome_annot'],
        # Convert list to string for the shell
        target_subtypes = ",".join(config['DG_Subset']['subset']['target_subtypes']),
        # Get the column name from config
        group_by = config['DG_Subset']['subset']['group_by'],
        bin_size = config['analysis_params']['tile_qc']['bin_size']
    output:
        subset = get_res_path(config['DG_Subset']['subset']['dir'] + "{sample}.h5ad")
    log: get_res_path(config['DG_Subset']['subset']['log'] + "{sample}.log")
    conda: '../envs/magic_env.yaml'
    shell:
        # Added {params.group_by} as the 5th argument
        "python scripts/subset_celltypes.py {input.adata} {wildcards.sample} {params.annotation} {params.target_subtypes} {params.group_by} {params.bin_size} {output.subset} > {log} 2>&1"

rule merge_DGsub:
    input: expand(get_res_path(config['DG_Subset']['subset']['dir'] + "{sample}.h5ad"), sample=SAMPLES)
    output:
        AnnDataSet = get_res_path(config['DG_Subset']['merge']['AnnDataSet']),
        flag = get_res_path(config['DG_Subset']['merge']['flag'])
    log: get_res_path(config['DG_Subset']['merge']['log'] + "merge.log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/3.Merge_AnnData.py {input} {output.AnnDataSet} {output.flag} > {log} 2>&1"

rule batch_DGsub:
    input: flag = get_res_path(config['DG_Subset']['merge']['flag'])
    params:
        dataset = get_res_path(config['DG_Subset']['merge']['AnnDataSet']),
        blacklist = config['path_to_blacklist'],
        n_feat = config['analysis_params']['batch']['n_features'],
        max_iter = config['analysis_params']['batch']['max_iter'],
        n_iter = config['analysis_params']['batch']['n_iter'],
        png = report(get_res_path(config['DG_Subset']['batch']['dir'] + "eigen.png"),category="CommonSamples",subcategory='BatchCorrection')
    output:
        flag_out = get_res_path(config['DG_Subset']['batch']['flag']),
        png1 = report(get_res_path(config['DG_Subset']['batch']['dir'] + "umap_sample.png"),category="CommonSamples",subcategory='BatchCorrection'),
        png2 = report(get_res_path(config['DG_Subset']['batch']['dir'] + "umap_aft_sample.png"),category="CommonSamples",subcategory='BatchCorrection'),
        png3 = report(get_res_path(config['DG_Subset']['batch']['dir'] + "umap_aft_leiden.png"),category="CommonSamples",subcategory='BatchCorrection')
    log: get_res_path(config['DG_Subset']['batch']['log'] + "batch.log")
    conda: '../envs/scanpy_env.yaml'
    shell:
        "python scripts/4.Batch_Correction.py {input.flag} {params.dataset} {params.blacklist} {params.n_feat} {params.max_iter} {params.n_iter} {params.png} {output.png1} {output.png2} {output.png3} {output.flag_out} > {log} 2>&1"

rule cluster_DGsub:
    input: flag = get_res_path(config['DG_Subset']['batch']['flag'])
    params:
        dataset = get_res_path(config['DG_Subset']['merge']['AnnDataSet']),
        res = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_res_path(config['DG_Subset']['clustering']['dir'])
    output:
        h5ad = get_res_path(config['DG_Subset']['clustering']['h5ad']),
        flag = get_res_path(config['DG_Subset']['clustering']['flag'])
    log: get_res_path(config['DG_Subset']['clustering']['log'] + "cluster.log")
    conda: '../envs/scanpy_env.yaml'
    shell:
        "python scripts/5.Clustering.py {input.flag} {params.dataset} {params.res} {params.out_dir} {output.h5ad} {output.flag} > {log} 2>&1"

rule gem_DGsub:
    input: h5ad = get_res_path(config['DG_Subset']['clustering']['h5ad'])
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
    output: gem = get_res_path(config['DG_Subset']['gem']['h5ad'])
    log: get_res_path(config['DG_Subset']['gem']['log'] + "gem.log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/6.Create_GEM.py {input.h5ad} {params.wd} {params.annot} {params.min_cells} {params.flavor} {params.batch} {params.n_top} {params.min_m} {params.max_m} {params.min_d} {params.n_jobs} {params.org} {params.rem_mit} {params.rem_ribo} {params.rem_sex} {params.coding} {params.meta} {output.gem} > {log} 2>&1
        """


rule create_peaks_DGsub:
    input: h5ad=get_res_path(config['DG_Subset']['clustering']['h5ad'])
    params:
        work_dir = config['DG_Subset']['peaks']['work_dir'],
        group_by = config['DG_Subset']['peaks']['group_by']
    output:
        peaks =  get_res_path(config['DG_Subset']['peaks']['h5ad'])
    log: get_res_path(config['DG_Subset']['peaks']['log'])
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/7.Create_Peaks.py {input.h5ad} {params.work_dir} {params.group_by} {output.peaks} > {log} 2>&1"