# rule targets:
#     input:
#         'results/{phenotype}/{group}.{genome}.report.md'

rule preprocess_data:
    input:
        R="code/preproc.R",
        csv="data/mikropml/{phenotype}/{group}.{genome}.csv"
    output:
        rds='data/mikropml/{phenotype}/{group}.{genome}.dat_proc.Rds'
    log:
        "log/{phenotype}/{group}.{genome}.preprocess_data.txt"
    benchmark:
        "benchmarks/{phenotype}/{group}.{genome}.preprocess_data.txt"
    params:
        outcome_colname='{phenotype}'
    resources:
        ncores=ncores, #This is present but not active in the cluster.json file since the multiple cores was not necessary for the size of input matrices, it was an inefficient use of resources; if preprocessing requires more, add back in this to the cluster.json file for this rule
        mem_mb = get_mem_mb_med
    script:
        "code/preproc.R"

rule run_ml:
    input:
        R="code/ml.R",
        rds=rules.preprocess_data.output.rds
    output:
        model="results/{phenotype}/runs/{group}.{genome}.{method}_{seed}_model.Rds",
        perf=temp("results/{phenotype}/runs/{group}.{genome}.{method}_{seed}_performance.csv"),
        feat=temp("results/{phenotype}/runs/{group}.{genome}.{method}_{seed}_feature-importance.csv")
    log:
        "log/runs/run_ml.{phenotype}.{group}.{genome}.{method}_{seed}.txt"
    benchmark:
        "benchmarks/runs/run_ml.{phenotype}.{group}.{genome}.{method}_{seed}.txt"
    params:
        outcome_colname='{phenotype}',
        method="{method}",
        seed="{seed}",
        kfold=kfold,
        alpha_var=alpha_var,
        lambda_var=lambda_var
    resources:
        ncores=ncores,
        mem_mb = get_mem_mb_low
    script:
        "code/ml.R"

rule combine_results:
    input:
        R="code/combine_results.R",
        csv=expand("results/{{phenotype}}/runs/{{group}}.{{genome}}.{method}_{seed}_{{type}}.csv", method=ml_methods, seed=seeds)
    output:
        csv='results/{phenotype}/{group}.{genome}.{type}_results.csv'
    log:
        "log/{phenotype}/{group}.{genome}.combine_results_{type}.txt"
    benchmark:
        "benchmarks/{phenotype}/{group}.{genome}.combine_results_{type}.txt"
    resources:
        mem_mb = get_mem_mb_lowest
    script:
        "code/combine_results.R"

rule combine_hp_performance:
    input:
        R='code/combine_hp_perf.R',
        rds=expand('results/{{phenotype}}/runs/{{group}}.{{genome}}.{{method}}_{seed}_model.Rds', seed=seeds)
    output:
        rds='results/{phenotype}/{group}.{genome}.hp_performance_results_{method}.Rds'
    log:
        "log/{phenotype}/{group}.{genome}.combine_hp_perf_{method}.txt"
    benchmark:
        "benchmarks/{phenotype}/{group}.{genome}.combine_hp_perf_{method}.txt"
    resources:
        mem_mb = get_mem_mb_med
    script:
        "code/combine_hp_perf.R"

rule combine_benchmarks:
    input:
        R='code/combine_benchmarks.R',
        # tsv=expand(rules.run_ml.benchmark, method = ml_methods, seed = seeds)
        tsv=expand("benchmarks/runs/run_ml.{{phenotype}}.{{group}}.{{genome}}.{method}_{seed}.txt", method=ml_methods, seed=seeds)
    output:
        csv='results/{phenotype}/{group}.{genome}.benchmarks_results.csv'
    resources:
        mem_mb = get_mem_mb_low
    log:
        'log/{phenotype}/{group}.{genome}.combine_benchmarks.txt'
    benchmark:
        "benchmarks/{phenotype}/{group}.{genome}.combine_benchmarks.txt"
    script:
        'code/combine_benchmarks.R'

rule plot_performance:
    input:
        R="code/plot_perf.R",
        csv='results/{phenotype}/{group}.{genome}.performance_results.csv'
    output:
        plot='figures/{phenotype}/{group}.{genome}.performance.png'
    resources:
        mem_mb = get_mem_mb_low
    log:
        "log/{phenotype}/{group}.{genome}.plot_performance.txt"
    script:
        "code/plot_perf.R"

rule plot_hp_performance:
    input:
        R='code/plot_hp_perf.R',
        rds=rules.combine_hp_performance.output.rds,
    output:
        plot='figures/{phenotype}/{group}.{genome}.hp_performance_{method}.png'
    resources:
        mem_mb = get_mem_mb_low
    log:
        'log/{phenotype}/{group}.{genome}.plot_hp_perf_{method}.txt'
    script:
        'code/plot_hp_perf.R'

rule plot_benchmarks:
    input:
        R='code/plot_benchmarks.R',
        csv=rules.combine_benchmarks.output.csv
    output:
        plot='figures/{phenotype}/{group}.{genome}.benchmarks.png'
    resources:
        mem_mb = get_mem_mb_low
    log:
        'log/{phenotype}/{group}.{genome}.plot_benchmarks.txt'
    script:
        'code/plot_benchmarks.R'

rule render_report:
    input:
        Rmd='report.Rmd',
        R='code/render.R',
        perf_plot=rules.plot_performance.output.plot,
        #hp_plot=expand(rules.plot_hp_performance.output.plot, method = ml_methods),
        hp_plot=expand('figures/{{phenotype}}/{{group}}.{{genome}}.hp_performance_{method}.png', method=ml_methods),
        bench_plot=rules.plot_benchmarks.output.plot
    output:
        doc='results/{phenotype}/{group}.{genome}.report.md'
    log:
        "log/{phenotype}/{group}.{genome}.render_report.txt"
    params:
        nseeds=nseeds,
        ml_methods=ml_methods,
        ncores=ncores,
        kfold=kfold,
    resources:
        mem_mb = "1GB"
    script:
        'code/render.R'

rule clean:
    input:
        rules.render_report.output,
        rules.plot_performance.output.plot,
        rules.plot_benchmarks.output.plot
    shell:
        '''
        rm -rf {input}
        '''
