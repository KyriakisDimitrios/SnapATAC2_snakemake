# --- 1. Standard Branch ---
rule merge_std:
    input:
        adatas = expand(
            get_path("standard", "feature_qc", "dir") + "{sample}.h5ad",
            sample=config['samples']
        ),
        annotation_gff3_file = config['genome_annot']
    output:
        AnnDataSet = get_path("standard", "merge", "AnnDataSet"),
        flag       = get_path("standard", "merge", "flag")
    log:
        get_path("standard", "merge", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/3.1.Merge_AnnData.py {input.annotation_gff3_file} {output.AnnDataSet} {output.flag}  {input.adatas} > {log} 2>&1"

# --- 2. Add Metadata (New Rule) ---
rule add_metadata:
    input:
        # Takes the output of the previous 'merge_std' rule
        merged_ds = get_path("standard", "merge", "AnnDataSet"),
        # Metadata files
        rush_meta = config['rush_metadata'],
        va_meta   = config['va_metadata']
    output:
        # Creates a NEW annotated file
        annotated_ds = get_path("standard", "add_metadata", "h5ad_output"),
    log:
        get_path("standard", "add_metadata", "log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/3.2.AddMetadata.py \
        {input.merged_ds} \
        {input.rush_meta} \
        {input.va_meta} \
        {output.annotated_ds} \
        > {log} 2>&1
        """

# --- 3. Filter Samples Branch (High Quality) ---
rule sample_qc_filter:
    input:
        # Input 1: The annotated single H5AD from Step 3.2
        annotated_ds=get_path("standard", "merge", "AnnDataSet"),#
        #annotated_ds=get_path("standard","add_metadata","h5ad_output"),

        # Input 2: The raw stats text files from Step 1
        # We use the path defined in 'structure' -> 'tile_qc' -> 'raw_qc'
        raw_stats=expand(
            get_path("standard","tile_qc","raw_qc"),
            sample=config['samples']
        )
    params:
        # Use the parameter from your config
        min_tsse=config['analysis_params']['sample_qc_filter']['min_mean_tsse']
    output:
        # Output: Single H5AD file (Merged_WM_RS.h5ad)
        filtered_dir = directory(get_path("standard","sample_qc_filter","out_dir")),
        filtered_ds=get_path("standard","sample_qc_filter","h5ad_output")
    log:
        get_path("standard","sample_qc_filter","log")
    conda: '../envs/magic_env.yaml'
    shell:
        """
        python scripts/3.3.Filter_Samples.py \
        {input.annotated_ds} \
        {output.filtered_dir} \
        {params.min_tsse} \
        {input.raw_stats} \
        > {log} 2>&1
        """





#
#
#
#
# # --- 2. Metadata Branch ---
# rule merge_meta:
#     input:
#         adatas = expand(
#             os.path.join(config['results_dir'], config['branches']['metadata']['constrain']['dir'], "{sample}.h5ad"),
#             sample=config['samples']
#         )
#     output:
#         AnnDataSet = get_path("metadata", "merge", "AnnDataSet"),
#         flag       = get_path("metadata", "merge", "flag")
#     log:
#         get_path("metadata", "merge", "log")
#     conda: '../envs/magic_env.yaml'
#     shell:
#         "python scripts/3.1.Merge_AnnData.py {input.adatas}  {output.AnnDataSet} {output.flag} > {log} 2>&1"
#
#
# # --- 3. DG_Subset Branch ---
# rule merge_DGsub:
#     input:
#         adatas = expand(
#             get_path("DG_Subset", "subset", "h5ad"),
#             sample=config['samples']
#         )
#     output:
#         AnnDataSet = get_path("DG_Subset", "merge", "AnnDataSet"),
#         flag       = get_path("DG_Subset", "merge", "flag")
#     log:
#         get_path("DG_Subset", "merge", "log")
#     conda: '../envs/magic_env.yaml'
#     shell:
#         "python scripts/3.1.Merge_AnnData.py {input.adatas} {output.AnnDataSet} {output.flag} > {log} 2>&1"