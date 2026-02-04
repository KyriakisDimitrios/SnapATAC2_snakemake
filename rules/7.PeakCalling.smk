# --- 1. Standard Branch Peak Calling ---
rule create_peaks_std:
    input:
        # Points to standard/05.clustering/Clustered.h5ad.gz
        h5ad = get_path("standard", "clustering", "h5ad")
    params:
        work_dir = get_path("standard", "peaks", "work_dir"),
        # Uses config['branches']['standard']['peaks']['group_by'] -> e.g., 'leiden'
        group_by = config['branches']['standard']['peaks']['group_by']
    output:
        peaks = get_path("standard", "peaks", "h5ad")
    log:
        get_path("standard", "peaks", "log")
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

# --- 2. Metadata Branch Peak Calling ---
rule create_peaks_meta:
    input:
        # Points to metadata/05.clustering/Clustered_Meta.h5ad.gz
        h5ad = get_path("metadata", "clustering", "h5ad")
    params:
        work_dir = get_path("metadata", "peaks", "work_dir"),
        # Uses config['branches']['metadata']['peaks']['group_by']
        group_by = config['branches']['metadata']['peaks']['group_by']
    output:
        peaks = get_path("metadata", "peaks", "h5ad")
    log:
        get_path("metadata", "peaks", "log")
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
