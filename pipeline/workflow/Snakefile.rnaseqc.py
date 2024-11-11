__author__ = "John Su"

#include: "path/to/other/snakefile" #relative to current file
#include: "rules/quant.smk"

configfile: workflow.basedir + "/../config/config.yaml"
scripts_folder=workflow.basedir + "/scripts"
resources_folder=workflow.basedir + "/../resources/"

#snakemake --cores all -F --use-singularity --singularity-args \"-B DIR_OUTSIDER:DIR_DOCKER\" -C out_dir=$PWD meta_dir=meta_ERP003613_SRP028336 -s Snakefile.XYZ.py TARGET_FILE

#extra = snakemake.params.get("extra", "")
#index = snakemake.input.idx
#config in command line: snakemake -npf --configfile=config.yaml

# please declare all variables related to external resources here, and specify them in config file or in -C arguement in command
#ref_genome=config["ref_genome"]
#reads1=config["reads1"]
#reads2=config["reads2"]
#Singuarlity image directory
#SI=config["SI"]

rnaseqc_gtf = resources_folder + "/gencode.v39.GRCh38.genes.gtf"


rule rnaseqc:
    input:
        Sorted_bam="STAR_Output/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "rnaseqc_out/{sample}/{sample}.Aligned.sortedByCoord.out.bam.metrics.tsv"
    log:
        "logs/rnaseqc/{sample}.log"
    params:
        ""
    threads:  # Samtools takes additional threads through its option -@
        2     # This value - 1 will be sent to -@
    shell:
        "mkdir -p rnaseqc_out/{wildcards.sample};rnaseqc " + rnaseqc_gtf + " {input.Sorted_bam} rnaseqc_out/{wildcards.sample}/ -u"
 
