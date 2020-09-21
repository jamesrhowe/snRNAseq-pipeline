SAMPLES = ["ASt-1", "ASt-2", "CeA-1", "CeA-2", "CeA-3", "CPb-1", "CPt-1", "CPt-2", "CPt-3"]
REGIONS = ["ASt", "CeA", "CPt"]


# rmarkdown::render('analysis/notebooks/Preprocessing.Rmd', params = list(input_path = "../../datasets/raw/Ast-1_raw", dataset = "Ast-1"))
# snakemake -p datasets/preprocessed/Ast-1_preprocessed.rds --cores 4
rule preprocessing:
    params:
        dataset='{sample}',
        input_path="../../datasets/raw/{sample}_raw/"
    output:
        "datasets/preprocessed/{sample}_preprocessed.rds"
    shell:
        """
        Rscript -e \"rmarkdown::render('analysis/notebooks/Preprocessing.Rmd', params = list(input_path = '{params.input_path}', dataset = '{params.dataset}'))\"
        """


# snakemake -p datasets/merged/ASt_preprocessed.rds datasets/merged/CeA_preprocessed.rds --cores 4
#rule merge:
#    output:
#        "datasets/merged/{smaple}_merged.rds"
#    params:
#        data="{sample}",
#        input="datasets/preprocessed/{regions}_raw/"
#    shell:
#       "Rscript --vanilla analysis/notebooks/MergeNormalize.R {params.input} {params.data}"
