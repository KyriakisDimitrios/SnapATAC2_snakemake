rule create_gem:
    input:
        h5ad_input = os.path.join(config['results_dir_path'], config['clustering']["output_h5ad_path"])
    params:
        work_dir = config["projdir"], # Using global projdir
        genome_annot = config["genome_annot"],
        min_cells = config["create_gem"]["min_cells"],
        flavor = config["create_gem"]["flavor"],
        batch_key = config["create_gem"]["batch_key"],
        n_top_genes = config["create_gem"]["n_top_genes"],
        min_mean = config["create_gem"]["min_mean"],
        max_mean = config["create_gem"]["max_mean"],
        min_disp = config["create_gem"]["min_disp"],
        n_jobs = config["create_gem"]["n_jobs"],
        organism = config["create_gem"]["organism"],
        remove_mit = config["create_gem"]["remove_mit"],
        remove_ribo = config["create_gem"]["remove_ribo"],
        remove_sex_genes = config["create_gem"]["remove_sex_genes"],
        only_coding_genes = config["create_gem"]["only_coding_genes"],
        metadata_file_path = config["ensemble_database"] # Fixed reference
    output:
        output_h5ad = os.path.join(config['results_dir_path'], config["create_gem"]["counts_output_h5ad"])
    log:
        log = os.path.join(config['results_dir_path'], config["create_gem"]["log"])
    shell:
        """
        source activate /sc/arion/projects/CommonMind/kyriad02/conda/envs/snapatac2Jan2025
        python workflow/scripts/6.Create_GEM.py \
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
            {output.output_h5ad} > {log} 2>&1
        """