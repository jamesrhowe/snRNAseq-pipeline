SAMPLES = ["ASt-1", "ASt-2", "CeA-1", "CeA-2", "CeA-3", "CPb-1", "CPt-1", "CPt-2", "CPt-3"]
REGIONS = ["ASt", "CeA", "CPt"]


# snakemake -p datasets/preprocessed/ASt-1_preprocessed.rds datasets/preprocessed/ASt-2_preprocessed.rds --cores 4
rule preprocessing:
    params:
        data="{sample}",
        input="datasets/raw/{sample}_raw/"
    output:
        "datasets/preprocessed/{sample}_preprocessed.rds"
    shell:
        "Rscript --vanilla analysis/notebooks/Preprocessing.R {params.input} {params.data}"


# snakemake -p datasets/merged/ASt_preprocessed.rds datasets/merged/CeA_preprocessed.rds --cores 4
#rule merge:
#    output:
#        "datasets/merged/{smaple}_merged.rds"
#    params:
#        data="{sample}",
#        input="datasets/preprocessed/{regions}_raw/"
#    shell:
#       "Rscript --vanilla analysis/notebooks/MergeNormalize.R {params.input} {params.data}"
