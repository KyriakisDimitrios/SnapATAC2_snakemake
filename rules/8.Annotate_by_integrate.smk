rule annotate_integration:
    input:
        query_h5ad=get_path("standard","gem","h5ad_output"),
        reference_h5ad=config['analysis_params']['annotate_by_integrate']['reference_h5ad']
    params:
        ref_celltype_col=config['analysis_params']['annotate_by_integrate']['ref_celltype_col'],
        min_cells=config['analysis_params']['annotate_by_integrate']['min_cells'],
        n_top_genes=config['analysis_params']['annotate_by_integrate']['n_top_genes'],
        n_layers=config['analysis_params']['annotate_by_integrate']['n_layers'],
        n_latent=config['analysis_params']['annotate_by_integrate']['n_latent'],
        max_epochs_scvi=config['analysis_params']['annotate_by_integrate']['max_epochs_scvi'],
        max_epochs_scanvi=config['analysis_params']['annotate_by_integrate']['max_epochs_scanvi'],
        n_samples_per_label=config['analysis_params']['annotate_by_integrate']['n_samples_per_label']
    output:
        csv_output=get_path("standard","annotate_by_integrate","csv_output"),
        flag_out=get_path("standard","annotate_by_integrate","flag"),

        # --- scVI / scANVI Diagnostics ---
        png_scvi_elbo=report(
            get_path("standard","annotate_by_integrate","dir") + "scvi_elbo.png",
            category="Annotation",
            subcategory="scVI Training"
        ),
        png_scanvi_elbo=report(
            get_path("standard","annotate_by_integrate","dir") + "scanvi_elbo.png",
            category="Annotation",
            subcategory="scANVI Training"
        ),
        png_umap_integration=report(
            get_path("standard","annotate_by_integrate","dir") + "umap_integration.png",
            category="Annotation",
            subcategory="Joint Embedding"
        )
    log:
        get_path("standard","annotate_by_integrate","log")
    conda: '../envs/scvi_integration_env.yaml'
    shell:
        """
        python scripts/8.Annotation_by_integrading.py \
        {input.reference_h5ad} \
        {input.query_h5ad} \
        {params.ref_celltype_col} \
        {params.min_cells} \
        {params.n_top_genes} \
        {params.n_layers} \
        {params.n_latent} \
        {params.max_epochs_scvi} \
        {params.max_epochs_scanvi} \
        {params.n_samples_per_label} \
        {output.png_scvi_elbo} \
        {output.png_scanvi_elbo} \
        {output.png_umap_integration} \
        {output.csv_output} \
        {output.flag_out} \
        > {log} 2>&1
        """