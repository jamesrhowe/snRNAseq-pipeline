configfile: "snakemake_config.yaml"

SAMPLE=config["sample"]
GROUP=config["group"]
THREADS=config["threads"]

rule all:
    input:
        expand("datasets/merged/{group}_merged.rds", group=GROUP)

#rule cluster:
#    input:
#        "datasets/merged/{group}_merged.rds"
#    conda:
#        "envs/rmarkdown.yaml"
#    log:
#        "logs/cluster_{group}.txt"
#    shell: """
#        Rscript -e "renv::restore()"
#        Rscript -e "rmarkdown::render(input = 'notebooks/Cluster_hicat.Rmd',
#                                      output_file = '../results/{wildcards.group}_cluster.nb.html',
#                                      params = list(title = 'Unsupervised Clustering: {wildcards.group}',
#                                                    group = '{wildcards.group}',
#                                                    input_path = '{input}',
#                                                    output_dir = '../datasets/merged/'))" 2> {log}
#    """

rule merge:
    input:
        "datasets/temp/{group}_list.rds"
    output:
        "datasets/merged/{group}_merged.rds",
        "results/{group}_merged.nb.html",
    conda:
        "envs/rmarkdown.yaml"
    log:
        "logs/merge_{group}.txt"
    shell: """
        Rscript -e "renv::restore()"
        Rscript -e "rmarkdown::render(input = 'notebooks/MergeNormalize.Rmd',
                                      output_file = '../results/{wildcards.group}_merged.nb.html',
                                      params = list(title = 'Merging and normalization: {wildcards.group}',
                                                    group = '{wildcards.group}',
                                                    input_path = '{input}',
                                                    output_dir = '../datasets/merged/'))" 2> {log}
        """

rule make_temp_merge_list:
    input:
        expand("datasets/preprocessed/{sample}_preprocessed.rds", sample=SAMPLE)
    output:
        temp("datasets/temp/{group}_list.rds")
    script:
        "scripts/make_merge_list.R"

rule preprocess:
    input:
        "datasets/raw/{sample}/matrix.mtx.gz"
    output:
        "datasets/preprocessed/{sample}_preprocessed.rds",
        "results/{sample}_preprocessed.nb.html"
    conda:
        "envs/rmarkdown.yaml"
    log:
        "logs/preprocess_{sample}.txt"
    shell: """
        Rscript -e "renv::restore()"
        Rscript -e "rmarkdown::render(input = 'notebooks/Preprocessing.Rmd',
                                      output_file = '../results/{wildcards.sample}_preprocessed.nb.html',
                                      params = list(title = 'Preprocessing: {wildcards.sample}',
                                                    dataset = '{wildcards.sample}',
                                                    input_path = '../datasets/raw/',
                                                    output_path = '../datasets/preprocessed/'))" 2> {log}
        """
