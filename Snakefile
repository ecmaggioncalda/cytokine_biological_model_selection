configfile: 'config/config.yml'

phenotype = config['phenotype']
genome = config['genome']

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']
alpha_var = config['alpha_var']
lambda_var = config['lambda_var']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

# attempt limit is set to 5 in the file code/submit_slurm.sh under restart-times
def get_mem_mb_lowest(wildcards, attempt):
    mem = attempt*1
    return "%dGB" % (mem)

def get_mem_mb_low(wildcards, attempt):
     mem = attempt*2
     return "%dGB" % (mem)

def get_mem_mb_med(wildcards, attempt):
     mem = attempt*4
     return "%dGB" % (mem)

def get_mem_mb_high(wildcards, attempt):
     mem = attempt*8
     return "%dGB" % (mem)

def get_geno_path(wildcards):
    return config["genome"][wildcards.genome]

rule all:
    input:
        expand("aggregated/{phenotype}.runs.csv", phenotype=phenotype)

# The group wildcard will define the different phenotype files generated during this checkpoint, which can be the single raw file, or multiple covariate files depending on user input in the prepro_overall.R script
# Once this checkpoint is complete, snakemake will re-evaluate the jobs that are required to complete the necessary downstream file creation
checkpoint prepro_overall:
    input:
        R = "code/prepro_overall.R",
        metadata = config['metadata']
    output:
        dat_dir = directory("data/pheno/{phenotype}")
    log:
        "log/{phenotype}.prepro_overall.txt"
    resources:
        ncores = ncores,
        mem_mb = get_mem_mb_lowest
    benchmark:
        "benchmarks/{phenotype}.prepro_overall.txt"
    script:
        "code/prepro_overall.R"

# This rule generates the phenotype:genotype data frame in the necessary input format for mikropml's preprocessing function (see rule in mikropml.smk, preprocess_data)
rule generate_mikropml_df:
    input:
        R = "code/generate_mikropml_df.R",
        pheno = "data/pheno/{phenotype}/{group}.tsv",
        feature_file = "data/minimal_filtered_features.csv"
    output:
        file_name = "data/mikropml/{phenotype}/{group}.{genome}.csv"
    params:
        path = get_geno_path
    log:
        "log/{phenotype}/{group}.{genome}.generate_mikropml_df.txt"
    resources:
        ncores = ncores,
        mem_mb = get_mem_mb_lowest
    benchmark:
        "benchmarks/{phenotype}/{group}.{genome}.generate_mikropml_df.txt"
    script:
        "code/generate_mikropml_df.R"

include: "mikropml.smk"

def aggregate_input4(wildcards):
    checkpoint_output = checkpoints.prepro_overall.get(**wildcards).output[0]
    return expand('results/{phenotype}/{group}.{genome}.report.md',
        phenotype=wildcards.phenotype,
        group=glob_wildcards(os.path.join(checkpoint_output,"{group}.tsv")).group,
        genome=genome)

# finish_list = [aggregate_input1, aggregate_input2, aggregate_input4]
# finish_list = [aggregate_input1, aggregate_input4]

rule finish_test:
    input:
        # finish_list
        aggregate_input4
    output:
        "aggregated/{phenotype}.runs.csv"
    log:
        "log/{phenotype}/finish.txt"
    resources:
        mem_mb = get_mem_mb_lowest
    script:
        "code/assemble_files.py"
