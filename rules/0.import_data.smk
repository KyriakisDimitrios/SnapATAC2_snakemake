rule import_data:
    input:
        frags = os.path.join(config['cellranger_dir'], "{sample}/outs/atac_fragments.tsv.gz"),
        metadata_path = config['metadata_path']
    params:
        min_num_fragments = config['import_data']['min_num_fragments'],
        n_jobs = config['import_data']['n_jobs']
    output:
        h5ad = os.path.join(config['results_dir_path'], config['import_data']['h5ads_output_dir'], "{sample}.h5ad")
    log:
        log = os.path.join(config['results_dir_path'], config['import_data']['log'], "0.import_data_{sample}.log")
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/0.import_data.py {input.frags} {input.metadata_path} {wildcards.sample} {params.min_num_fragments} {params.n_jobs} {output.h5ad}  > {log} 2>&1"