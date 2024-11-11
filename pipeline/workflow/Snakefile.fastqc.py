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
ref_genome=config["ref_genome"]
#reads1=config["reads1"]
#reads2=config["reads2"]
#Singuarlity image directory
SI=config["SI"]



rule fastqc_pe:
    input:
        reads1="fq/{sample}_1.fq.gz",
        reads2="fq/{sample}_2.fq.gz"
    output:
        "fastqc_out/{sample}/{sample}_2_fastqc.html"
    log:
        "logs/fastqc_out/{sample}.log"
    params:
        ""
    threads:  # Samtools takes additional threads through its option -@
        2     # This value - 1 will be sent to -@
    shell:
        "mkdir -p fastqc_out/{wildcards.sample} $TMPDIR/fastqc;fastqc -o fastqc_out/{wildcards.sample}/ --extract -f fastq -d $TMPDIR/fastqc {input.reads1} {input.reads2}"

rule fastqc_se:
    input:
        reads1="fq/{sample}_1.fq.gz"
    output:
        "fastqc_out/{sample}/{sample}_1_fastqc.html"
    log:
        "logs/fastqc_out/{sample}.log"
    params:
        ""
    threads:  # Samtools takes additional threads through its option -@
        2     # This value - 1 will be sent to -@
    shell:
        "mkdir -p fastqc_out/{wildcards.sample} $TMPDIR/fastqc;fastqc -o fastqc_out/{wildcards.sample}/ --extract -f fastq -d $TMPDIR/fastqc {input.reads1}" 
