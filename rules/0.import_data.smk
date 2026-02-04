# Function to determine which fragment file exists
def get_fragment_file(wildcards):
    # Construct both possible paths
    base_path = os.path.join(config['cellranger_dir'],wildcards.sample,"outs")

    option_1 = os.path.join(base_path,"fragments.tsv.gz")
    option_2 = os.path.join(base_path,"atac_fragments.tsv.gz")

    # Check which one exists
    if os.path.exists(option_2):
        return option_2

    # Default to option_1 (if neither exists, Snakemake will show a "Missing Input" error for this path)
    return option_1


rule import_data:
    input:
        frags=get_fragment_file  # <--- Pass the function name, don't call it ()
    params:
        min_num_fragments=config['import_data']['min_num_fragments'],
        n_jobs=config['import_data']['n_jobs']
    output:
        h5ad=os.path.join(config['results_dir_path'],config['import_data']['h5ads_output_dir'],"{sample}.h5ad")
    log:
        log=os.path.join(config['results_dir_path'],config['import_data']['log'])
    conda: '../envs/magic_env.yaml'
    shell:
        "python scripts/0.import_data.py {input.frags} {wildcards.sample} {params.min_num_fragments} {params.n_jobs} {output.h5ad} > {log} 2>&1"



# rule import_data:
#     input:
#         frags = os.path.join(config['cellranger_dir'], "{sample}/outs/fragments.tsv.gz")
#     params:
#         min_num_fragments = config['import_data']['min_num_fragments'],
#         n_jobs = config['import_data']['n_jobs']
#     output:
#         h5ad = os.path.join(config['results_dir_path'], config['import_data']['h5ads_output_dir'], "{sample}.h5ad")
#     log:
#         log = os.path.join(config['results_dir_path'], config['import_data']['log'], "0.import_data_{sample}.log")
#     conda: '../envs/magic_env.yaml'
#     shell:
#         "python scripts/0.import_data.py {input.frags} {wildcards.sample} {params.min_num_fragments} {params.n_jobs} {output.h5ad}  > {log} 2>&1"