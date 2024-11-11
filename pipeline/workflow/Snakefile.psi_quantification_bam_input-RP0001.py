__author__ = "John Su"

#import os
include: "Snakefile.rnaseqc.py"
include: "Snakefile.fastqc.py"
configfile: workflow.basedir + "/../config/config.yaml"
scripts_folder=workflow.basedir + "/scripts"
resources_folder=workflow.basedir + "/../resources/"

#snakemake --cores all -F --use-singularity --singularity-args \"-B DIR_OUTSIDER:DIR_DOCKER\" -C out_dir=$PWD meta_dir=meta_ERP003613_SRP028336 -s Snakefile.XYZ.py TARGET_FILE

#extra = snakemake.params.get("extra", "")
#index = snakemake.input.idx
#config in command line: snakemake -npf --configfile=config.yaml
#or specify this in this file: configfile: "config.yaml"
ref_genome=config["ref_genome"]
bam=config["bam"]
read_length=config["read_length"]
#Singuarlity image directory
SI=config["SI"]



rule bam2fastq:
    input:
        bam=bam
    output:
        reads1=temp("fq/{sample}_1.fq.gz"),
        reads2=temp("fq/{sample}_2.fq.gz"),
        nsorted_bam=temp("fq/{sample}_nsorted.bam"),
        se_reads=temp("fq/{sample}_SE.fq.gz")
    log:
        "logs/fq/{sample}.log"
    params:
        ""# optional params string
    threads:  # Samtools takes additional threads through its option -@
        10     # This value - 1 will be sent to -@
    shell:
        "samtools sort -n {input.bam} -o {output.nsorted_bam}; samtools fastq -1 {output.reads1} -2 {output.reads2} -s {output.se_reads} {output.nsorted_bam}"

rule align:
    input:
        reads1="fq/{sample}_1.fq.gz",
        reads2="fq/{sample}_2.fq.gz"
    output:
        SJ_tab=temp("STAR_Output/{sample}.SJ.out.tab"),
        Aligned_bam=temp("STAR_Output/{sample}.Aligned.out.bam"),
        Sorted_bam=temp("STAR_Output/{sample}.Aligned.sortedByCoord.out.bam")
    log:
        "logs/align/{sample}.log"
    params:
        ""# optional params string
    threads:  # Samtools takes additional threads through its option -@
        10     # This value - 1 will be sent to -@
    shell:
        "rm -rf $TMPDIR/star;STAR --runThreadN 10  --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted SortedByCoordinate --genomeDir " + ref_genome + " --outTmpDir $TMPDIR/star --readFilesIn {input.reads1} {input.reads2} --outFileNamePrefix STAR_Output/{wildcards.sample}. --readFilesCommand zcat" 
        #"rm -rf $TMPDIR/star;STAR --runThreadN 10  --outSAMtype BAM Unsorted SortedByCoordinate --genomeDir " + ref_genome + " --outTmpDir $TMPDIR/star --readFilesIn {input.reads1} {input.reads2} --outFileNamePrefix STAR_Output/{wildcards.sample}. --readFilesCommand zcat" 

rule cal_psi:
    input:
        SJ_tab="STAR_Output/{sample}.SJ.out.tab",
        Aligned_bam="STAR_Output/{sample}.Aligned.out.bam"
    output:
        junction_bed=temp("PSI_Calculation/{sample}.junctions.bed"),
        incl=temp("PSI_Calculation/{sample}.PSI_exonic_parts.inclusion"),
        excl=temp("PSI_Calculation/{sample}.PSI_exonic_parts.exclusion"),
        gff=temp("PSI_Calculation/{sample}.PSI_exonic_parts.gff"),
        PSI="PSI_Calculation/{sample}.PSI_exonic_parts.psi"

    log:
        "logs/cal_psi/{sample}.log"
    params:
        ""
    container:
        SI + "/jvivian_bedtools.sif"# optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    shell:
        scripts_folder + "/generate_junction_bed.sh {input.SJ_tab} {output.junction_bed};" + scripts_folder + "/PSI.sh StartPSI " + resources_folder + "/gencode.v38.annotation.Exonic_Parts.LRG.noDup.gtf " + str(read_length) + " {input.Aligned_bam} {output.junction_bed} PSI_Calculation/{wildcards.sample}.PSI"
        
rule do_qc_and_psi:
    input:
        "PSI_Calculation/{sample}.PSI_exonic_parts.psi",
        "rnaseqc_out/{sample}/{sample}.Aligned.sortedByCoord.out.bam.metrics.tsv",
        "fastqc_out/{sample}/{sample}_2_fastqc.html"
    output:
        "PSI_Calculation/{sample}.DONE"
    shell:
        "touch PSI_Calculation/{wildcards.sample}.DONE"

rule do_qc_only:
    input:
        "rnaseqc_out/{sample}/{sample}.Aligned.sortedByCoord.out.bam.metrics.tsv",
        "fastqc_out/{sample}/{sample}_2_fastqc.html"
    output:
        "PSI_Calculation/{sample}.QC"
    shell:
        "touch PSI_Calculation/{wildcards.sample}.QC"


