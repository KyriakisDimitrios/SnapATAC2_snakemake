# --- 1. Standard Branch ---
rule gem_std:
    input:
        h5ad_input = get_path("standard","batch","h5ad_output") #get_path("standard", "clustering", "h5ad")
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
        output_h5ad = get_path("standard", "gem", "h5ad")
    log:
        get_path("standard", "gem", "log")
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
#
# # --- 2. Metadata Branch ---
# rule gem_meta:
#     input:
#         h5ad_input = get_path("metadata", "clustering", "h5ad")
#     params:
#         work_dir      = config["projdir"],
#         genome_annot  = config["genome_annot"],
#         min_cells     = config["analysis_params"]["gem"]["min_cells"],
#         flavor        = config["analysis_params"]["gem"]["flavor"],
#         batch_key     = config["analysis_params"]["gem"]["batch_key"],
#         n_top_genes   = config["analysis_params"]["gem"]["n_top_genes"],
#         min_mean      = config["analysis_params"]["gem"]["min_mean"],
#         max_mean      = config["analysis_params"]["gem"]["max_mean"],
#         min_disp      = config["analysis_params"]["gem"]["min_disp"],
#         n_jobs        = config["analysis_params"]["gem"]["n_jobs"],
#         organism      = config["analysis_params"]["gem"]["organism"],
#         remove_mit       = config["analysis_params"]["gem"]["remove_mit"],
#         remove_ribo      = config["analysis_params"]["gem"]["remove_ribo"],
#         remove_sex_genes = config["analysis_params"]["gem"]["remove_sex_genes"],
#         only_coding_genes= config["analysis_params"]["gem"]["only_coding_genes"],
#         metadata_file_path = config["ensemble_database"]
#     output:
#         output_h5ad = get_path("metadata", "gem", "h5ad")
#     log:
#         get_path("metadata", "gem", "log")
#     conda: '../envs/magic_env.yaml'
#     shell:
#         """
#         python scripts/6.Create_GEM.py \
#             {input.h5ad_input} \
#             {params.work_dir} \
#             {params.genome_annot} \
#             {params.min_cells} \
#             {params.flavor} \
#             {params.batch_key} \
#             {params.n_top_genes} \
#             {params.min_mean} \
#             {params.max_mean} \
#             {params.min_disp} \
#             {params.n_jobs} \
#             {params.organism} \
#             {params.remove_mit} \
#             {params.remove_ribo} \
#             {params.remove_sex_genes} \
#             {params.only_coding_genes} \
#             {params.metadata_file_path} \
#             {output.output_h5ad} \
#             > {log} 2>&1
#         """
#
# # --- 3. DG Subset Branch ---
# rule gem_DGsub:
#     input:
#         h5ad_input = get_path("DG_Subset", "clustering", "h5ad")
#     params:
#         work_dir      = config["projdir"],
#         genome_annot  = config["genome_annot"],
#         min_cells     = config["analysis_params"]["gem"]["min_cells"],
#         flavor        = config["analysis_params"]["gem"]["flavor"],
#         batch_key     = config["analysis_params"]["gem"]["batch_key"],
#         n_top_genes   = config["analysis_params"]["gem"]["n_top_genes"],
#         min_mean      = config["analysis_params"]["gem"]["min_mean"],
#         max_mean      = config["analysis_params"]["gem"]["max_mean"],
#         min_disp      = config["analysis_params"]["gem"]["min_disp"],
#         n_jobs        = config["analysis_params"]["gem"]["n_jobs"],
#         organism      = config["analysis_params"]["gem"]["organism"],
#         remove_mit       = config["analysis_params"]["gem"]["remove_mit"],
#         remove_ribo      = config["analysis_params"]["gem"]["remove_ribo"],
#         remove_sex_genes = config["analysis_params"]["gem"]["remove_sex_genes"],
#         only_coding_genes= config["analysis_params"]["gem"]["only_coding_genes"],
#         metadata_file_path = config["ensemble_database"]
#     output:
#         output_h5ad = get_path("DG_Subset", "gem", "h5ad")
#     log:
#         get_path("DG_Subset", "gem", "log")
#     conda: '../envs/magic_env.yaml'
#     shell:
#         """
#         python scripts/6.Create_GEM.py \
#             {input.h5ad_input} \
#             {params.work_dir} \
#             {params.genome_annot} \
#             {params.min_cells} \
#             {params.flavor} \
#             {params.batch_key} \
#             {params.n_top_genes} \
#             {params.min_mean} \
#             {params.max_mean} \
#             {params.min_disp} \
#             {params.n_jobs} \
#             {params.organism} \
#             {params.remove_mit} \
#             {params.remove_ribo} \
#             {params.remove_sex_genes} \
#             {params.only_coding_genes} \
#             {params.metadata_file_path} \
#             {output.output_h5ad} \
#             > {log} 2>&1
#         """