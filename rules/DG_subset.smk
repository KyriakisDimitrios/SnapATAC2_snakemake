# ==============================================================================
# BRANCH B: METADATA CONSTRAINED (Whitelist)
# Flow: Import -> Constrain -> Merge -> Batch -> Cluster -> GEM
# ==============================================================================

rule subset_DGsub:
    input:
        # Input: 00.raw_fragments/{sample}.h5ad
        # We join results path with the import_data relative path
        adata = os.path.join(config['results_dir_path'], config['import_data']['h5ads_output'])
    params:
        annotation = config['genome_annot'],
        # Join list into string for shell argument
        target_subtypes = ",".join(config['branches']['DG_Subset']['subset']['target_subtypes']),
        group_by = config['branches']['DG_Subset']['subset']['group_by'],
        bin_size = config['analysis_params']['tile_qc']['bin_size']
    output:
        # FIX: Use 'h5ad_out' which contains "{sample}_subsetDG.h5ad" defined in config
        output_file = get_path("DG_Subset", "subset_DGsub", "h5ad_out")
    log:
        get_path("DG_Subset", "subset_DGsub", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/subset_celltypes.py \
        {input.adata} \
        {wildcards.sample} \
        {params.annotation} \
        {params.target_subtypes} \
        {params.group_by} \
        {params.bin_size} \
        {output.output_file} \
        > {log} 2>&1
        """

rule merge_DGsub:
    input:
        adatas = expand(get_path("DG_Subset", "subset_DGsub", "h5ad_out"),sample=config['samples']),
        metadata_path= config['metadata_path']
    output:
        AnnDataSet = get_path("DG_Subset","merge","AnnDataSet"),
    log: get_path("DG_Subset", "merge", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/3.Merge_AnnData.py {input.metadata_path} {output.AnnDataSet}  {input.adatas}  > {log} 2>&1"


# --- 1. Standard Branch ---
rule batch_DGsub:
    input:
        h5ad_input   =  get_path("DG_Subset","merge","AnnDataSet")
    params:
        batch_var = config['analysis_params']['batch']['batch_var'],
        blacklist = config['path_to_blacklist'],
        n_feat    = config['analysis_params']['batch']['n_features'],
        max_iter  = config['analysis_params']['batch']['max_iter'],
        n_iter    = config['analysis_params']['batch']['n_iter'],
        # Eigen Plot (Param because it is an output of the script but not a rule output used downstream)
        png_eigen = report(get_path("DG_Subset", "batch", "dir") + "eigen.png", category="Standard", subcategory="Batch")
    output:
        h5ad_output = get_path("DG_Subset","batch","h5ad_output"),
        flag_out = get_path("DG_Subset", "batch", "flag"),
        png1     = report(get_path("DG_Subset", "batch", "dir") + "umap_sample.png", category="Standard", subcategory="Batch"),
        png2     = report(get_path("DG_Subset", "batch", "dir") + "umap_aft_sample.png", category="Standard", subcategory="Batch"),
        png3     = report(get_path("DG_Subset", "batch", "dir") + "umap_aft_leiden.png", category="Standard", subcategory="Batch")
    log:
        get_path("DG_Subset", "batch", "log")
    conda: '../envs/scanpy_env.yaml' # Assuming magic_env contains snapatac2/scanpy
    shell:
        """
        python scripts/4.Batch_Correction.py \
        {input.h5ad_input} \
        {params.batch_var} \
        {params.blacklist} \
        {params.n_feat} \
        {params.max_iter} \
        {params.n_iter} \
        {params.png_eigen} \
        {output.png1} \
        {output.png2} \
        {output.png3} \
        {output.flag_out} \
        {output.h5ad_output} \
        > {log} 2>&1
        """

# --- 3. DG Subset Branch ---
rule cluster_DGsub:
    input:
        flag = get_path("DG_Subset", "batch", "flag")
    params:
        dataset = get_path("DG_Subset", "merge", "AnnDataSet"),
        res     = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_path("DG_Subset", "clustering", "dir")
    output:
        h5ad = get_path("DG_Subset", "clustering", "h5ad"),
        flag = get_path("DG_Subset", "clustering", "flag")
    log:
        get_path("DG_Subset", "clustering", "log")
    conda: '../envs/scanpy_env.yaml'
    shell:
        """
        python scripts/5.Clustering.py \
        {input.flag} \
        {params.dataset} \
        {params.res} \
        {params.out_dir} \
        {output.h5ad} \
        {output.flag} \
        > {log} 2>&1
        """

# --- 3. DG Subset Branch ---
rule gem_DGsub:
    input:
        h5ad_input = get_path("DG_Subset","clustering","h5ad")
    params:
        work_dir=config["projdir"],
        genome_annot=config["genome_annot"],
        min_cells=config["analysis_params"]["gem"]["min_cells"],
        flavor=config["analysis_params"]["gem"]["flavor"],
        batch_key=config["analysis_params"]["gem"]["batch_key"],
        n_top_genes=config["analysis_params"]["gem"]["n_top_genes"],
        min_mean=config["analysis_params"]["gem"]["min_mean"],
        max_mean=config["analysis_params"]["gem"]["max_mean"],
        min_disp=config["analysis_params"]["gem"]["min_disp"],
        n_jobs=config["analysis_params"]["gem"]["n_jobs"],
        organism=config["analysis_params"]["gem"]["organism"],
        remove_mit=config["analysis_params"]["gem"]["remove_mit"],
        remove_ribo=config["analysis_params"]["gem"]["remove_ribo"],
        remove_sex_genes=config["analysis_params"]["gem"]["remove_sex_genes"],
        only_coding_genes=config["analysis_params"]["gem"]["only_coding_genes"],
        metadata_file_path=config["ensemble_database"]
    output:
        output_h5ad=get_path("DG_Subset","gem","h5ad")
    log:
        get_path("DG_Subset","gem","log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/6.Create_GEM.py \
            {input.h5ad_input} \
            {params.work_dir} \
            {params.genome_annot} \
            {params.min_cells} \
            {params.flavor} \
            {params.batch_key} \
            {params.n_top_genes} \
            {params.min_mean} \
            {params.max_mean} \
            {params.min_disp} \
            {params.n_jobs} \
            {params.organism} \
            {params.remove_mit} \
            {params.remove_ribo} \
            {params.remove_sex_genes} \
            {params.only_coding_genes} \
            {params.metadata_file_path} \
            {output.output_h5ad} \
            > {log} 2>&1
        """


# --- 3. DG_Subset Branch Peak Calling ---
rule create_peaks_DGsub:
    input:
        # Input is the Clustered Object (contains fragments & metadata)
        h5ad = get_path("DG_Subset", "clustering", "h5ad")
    params:
        # Path to the work directory (for MACS3 temp files)
        work_dir = get_path("DG_Subset", "peaks", "work_dir"),
        # The column to group by (defined in branches -> DG_Subset -> peaks)
        group_by = config['branches']['DG_Subset']['peaks']['group_by']
    output:
        peaks = get_path("DG_Subset", "peaks", "h5ad")
    log:
        get_path("DG_Subset", "peaks", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/7.Create_Peaks.py \
        {input.h5ad} \
        {params.work_dir} \
        {params.group_by} \
        {output.peaks} \
        > {log} 2>&1
        """

# rule cluster_DGsub:
#     input: flag = get_res_path(config['DG_Subset']['batch']['flag'])
#     params:
#         dataset = get_res_path(config['DG_Subset']['merge']['AnnDataSet']),
#         res = config['analysis_params']['clustering']['resolutions'],
#         out_dir = get_res_path(config['DG_Subset']['clustering']['dir'])
#     output:
#         h5ad = get_res_path(config['DG_Subset']['clustering']['h5ad']),
#         flag = get_res_path(config['DG_Subset']['clustering']['flag'])
#     log: get_res_path(config['DG_Subset']['clustering']['log'])
#     conda: '../envs/scanpy_env.yaml'
#     shell:
#         "python scripts/5.Clustering.py {input.flag} {params.dataset} {params.res} {params.out_dir} {output.h5ad} {output.flag} > {log} 2>&1"

# rule gem_DGsub:
#     input: h5ad = get_res_path(config['DG_Subset']['clustering']['h5ad'])
#     params:
#         wd = config['projdir'],
#         annot = config['genome_annot'],
#         min_cells = config['analysis_params']['gem']['min_cells'],
#         flavor = config['analysis_params']['gem']['flavor'],
#         batch = config['analysis_params']['gem']['batch_key'],
#         n_top = config['analysis_params']['gem']['n_top_genes'],
#         min_m = config['analysis_params']['gem']['min_mean'],
#         max_m = config['analysis_params']['gem']['max_mean'],
#         min_d = config['analysis_params']['gem']['min_disp'],
#         n_jobs = config['analysis_params']['gem']['n_jobs'],
#         org = config['analysis_params']['gem']['organism'],
#         rem_mit = config['analysis_params']['gem']['remove_mit'],
#         rem_ribo = config['analysis_params']['gem']['remove_ribo'],
#         rem_sex = config['analysis_params']['gem']['remove_sex_genes'],
#         coding = config['analysis_params']['gem']['only_coding_genes'],
#         meta = config['ensemble_database']
#     output: gem = get_res_path(config['DG_Subset']['gem']['h5ad'])
#     log: get_res_path(config['DG_Subset']['gem']['log'])
#     conda: '../envs/magic_env.yaml'
#     shell:
#         """
#         python scripts/6.Create_GEM.py {input.h5ad} {params.wd} {params.annot} {params.min_cells} {params.flavor} {params.batch} {params.n_top} {params.min_m} {params.max_m} {params.min_d} {params.n_jobs} {params.org} {params.rem_mit} {params.rem_ribo} {params.rem_sex} {params.coding} {params.meta} {output.gem} > {log} 2>&1
#         """
#
#
# rule create_peaks_DGsub:
#     input: h5ad=get_res_path(config['DG_Subset']['clustering']['h5ad'])
#     params:
#         work_dir = config['DG_Subset']['peaks']['work_dir'],
#         group_by = config['DG_Subset']['peaks']['group_by']
#     output:
#         peaks =  get_res_path(config['DG_Subset']['peaks']['h5ad'])
#     log: get_res_path(config['DG_Subset']['peaks']['log'])
#     conda: '../envs/magic_env.yaml'
#     shell:
#         "python scripts/7.Create_Peaks.py {input.h5ad} {params.work_dir} {params.group_by} {output.peaks} > {log} 2>&1"