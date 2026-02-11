branch_to_run = 'Subset_byRNA'


rule constrain_meta:
    input:
        adata=os.path.join(config['results_dir_path'],config['import_data']['h5ads_output'])
    params:
        annotation = config['genome_annot'],
        bin_size = config['analysis_params']['tile_qc']['bin_size']
    output:
        output_file = get_path(branch_to_run, "subset_byRNA", "h5ad_out")
    log: get_path(branch_to_run, "subset_byRNA", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/constrain_data.py {input.adata} {wildcards.sample} {params.annotation} {params.bin_size} {output.output_file} > {log} 2>&1"

rule merge_CTbyRNA:
    input:
        adatas = expand(get_path(branch_to_run, "subset_byRNA", "h5ad_out"),sample=config['samples']),
        metadata_path= config['metadata_path']
    output:
        AnnDataSet = get_path(branch_to_run,"merge","AnnDataSet"),
    log: get_path(branch_to_run, "merge", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/3.Merge_AnnData.py {input.metadata_path} {output.AnnDataSet}  {input.adatas}  > {log} 2>&1"


# --- 1. Standard Branch ---
rule batch_CTbyRNA:
    input:
        h5ad_input = get_path(branch_to_run,"merge","AnnDataSet")
    params:
        batch_var=config['analysis_params']['batch']['batch_var'],
        blacklist=config['path_to_blacklist'],
        n_feat=config['analysis_params']['batch']['n_features'],
        max_iter=config['analysis_params']['batch']['max_iter'],
        n_iter=config['analysis_params']['batch']['n_iter'],

        # --- 1. Diagnostics ---
        png_eigen=report(
            get_path(branch_to_run,"batch","dir") + "eigen.png",
            category="Batch Correction",
            subcategory="Diagnostics"
        ),

        # --- 2. Pre-Correction ---
        png_pre=report(
            get_path(branch_to_run,"batch","dir") + "umap_pre_sample.png",
            category="Batch Correction",
            subcategory="Pre-Correction"
        ),

        # --- 3. MNC Method ---
        png_mnc_s=report(
            get_path(branch_to_run,"batch","dir") + "umap_mnc_sample.png",
            category="Batch Correction",
            subcategory="MNC"
        ),
        png_mnc_l=report(
            get_path(branch_to_run,"batch","dir") + "umap_mnc_leiden.png",
            category="Batch Correction",
            subcategory="MNC"
        ),

        # --- 4. Harmony Method ---
        png_harm_s=report(
            get_path(branch_to_run,"batch","dir") + "umap_harmony_sample.png",
            category="Batch Correction",
            subcategory="Harmony"
        ),
        png_harm_l=report(
            get_path(branch_to_run,"batch","dir") + "umap_harmony_leiden.png",
            category="Batch Correction",
            subcategory="Harmony"
        )

    output:
        flag_out=get_path(branch_to_run,"batch","flag"),
        h5ad_output=get_path(branch_to_run,"batch","h5ad_output")
    log:
        get_path(branch_to_run,"batch","log")
    conda: '../envs/snapatac2_env.yaml'
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
        {params.png_pre} \
        {params.png_mnc_s} \
        {params.png_mnc_l} \
        {params.png_harm_s} \
        {params.png_harm_l} \
        {output.flag_out} \
        {output.h5ad_output} \
        > {log} 2>&1
        """

# --- 3. DG Subset Branch ---
rule cluster_CTbyRNA:
    input:
        h5ad_input = get_path(branch_to_run,"batch","h5ad_output"),
    params:
        res     = config['analysis_params']['clustering']['resolutions'],
        out_dir = get_path(branch_to_run, "clustering", "dir")
    output:
        h5ad_output = get_path(branch_to_run, "clustering", "h5ad"),
    log:
        get_path(branch_to_run, "clustering", "log")
    conda: '../envs/scanpy_env.yaml'
    shell:
        """
        python scripts/5.Clustering.py \
        {input.h5ad_input} \
        {params.res} \
        {params.out_dir} \
        {output.h5ad_output} > {log} 2>&1
        """


# --- 1. Standard Branch ---
rule gem_CTbyRNA:
    input:
        h5ad_input = get_path(branch_to_run,"batch","h5ad_output")
    params:
        work_dir      = config["projdir"],
        genome_annot  = config["genome_annot"],
        # Analysis Params
        min_cells     = config["analysis_params"]["gem"]["min_cells"],
        flavor        = config["analysis_params"]["gem"]["flavor"],
        batch_key     = config["analysis_params"]["gem"]["batch_key"],
        n_top_genes   = config["analysis_params"]["gem"]["n_top_genes"],
        min_mean      = config["analysis_params"]["gem"]["min_mean"],
        max_mean      = config["analysis_params"]["gem"]["max_mean"],
        min_disp      = config["analysis_params"]["gem"]["min_disp"],
        n_jobs        = config["analysis_params"]["gem"]["n_jobs"],
        organism      = config["analysis_params"]["gem"]["organism"],
        # Booleans (True/False)
        remove_mit       = config["analysis_params"]["gem"]["remove_mit"],
        remove_ribo      = config["analysis_params"]["gem"]["remove_ribo"],
        remove_sex_genes = config["analysis_params"]["gem"]["remove_sex_genes"],
        only_coding_genes= config["analysis_params"]["gem"]["only_coding_genes"],
        # Reference Database
        metadata_file_path = config["ensemble_database"]
    output:
        output_h5ad = get_path(branch_to_run, "gem", "h5ad")
    log:
        get_path(branch_to_run, "gem", "log")
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
rule create_peaks_CTbyRNA:
    input:
        # Input is the Clustered Object (contains fragments & metadata)
        h5ad = get_path(branch_to_run,"batch","h5ad_output")
    params:
        # Path to the work directory (for MACS3 temp files)
        work_dir = get_path(branch_to_run, "peaks", "work_dir"),
        # The column to group by (defined in branches -> DG_Subset -> peaks)
        group_by = config['branches'][branch_to_run]['peaks']['group_by']
    output:
        peaks = get_path(branch_to_run, "peaks", "h5ad")
    log:
        get_path(branch_to_run, "peaks", "log")
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
