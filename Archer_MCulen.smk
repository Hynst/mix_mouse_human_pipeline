# Archer panel - Martin Culen vzorky mix (clovek a mys)
#
import re
import datetime

import boto3
from snakemake.utils import min_version
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

configfile: ""
#container: "cerit.io/snakemake:v6.10.0"

# check minimal version
min_version("5.18.0")

# S3 credentials
AWS_ID=""
AWS_KEY=""
S3_BUCKET =""

# conect to remote S3
S3 = S3RemoteProvider(host="", access_key_id=AWS_ID, secret_access_key=AWS_KEY)

TARGET_BED = S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.bed")
TARGET_INTERVALS = TARGET_INTERVALS = S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.intervals")

# tools
CUTADAPT = "cutadapt"
BWA = "bwa"
SAMBLASTER= "samblaster"
SAMTOOLS= "samtools"
SAMBAMBA= "sambamba"
PICARD = S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/picard.jar")
UMI_TOOLS = "umi_tools"
MULTIQC = "multiqc"

GATK2= S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/GATK/gatk")

VEP= "vep"
BEDTOOLS= "bedtools"
VCF2TABLE= S3.remote(S3_BUCKET + "/src/myelo_proliferation/archer/reGenotype/vcf2table_MDS.py")
BCFTOOLS = "bcftools"

def vcfvep_list(wildcards):
        vcfvep_list = []
        run = config["Run"]
        samples = config["Samples"]
        for sample in samples:
            vcfvep_list = vcfvep_list + S3.remote([S3_BUCKET + "/data/mculen/" + run + "/samples/" + sample + "/variants/mutect2/" + sample \
            + ".mutect2.filt.norm.vep.csv"])
        return vcfvep_list

def covlist(wildcards):
    covlist = []
    run = config["Run"]
    samples = config["Samples"]
    for sample in samples:
            covlist = covlist + S3.remote([S3_BUCKET + "/data/mculen/" + run + "/samples/" + sample + "/coverage/" + sample \
            + ".perexon_stat.txt"])
    return covlist

def getFasta(wildcards):
    files = {'fwd': S3.remote(S3_BUCKET + config["Samples"][wildcards.sample].get('fwd'))}
    files['rev'] = S3.remote(S3_BUCKET + config["Samples"][wildcards.sample].get('rev'))
    return files

rule all:
    output:
        "../run_" + str(datetime.datetime.now()).replace(" ","-")
    input:
        covlist,
        vcfvep_list
    shell:
        """
        touch {output}
        """

########################################
#### VARIANT CALLING AND ANNOTATION ####
########################################
rule vcf2csv:
    output:
        csv = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vep.csv")
    input:
        vcf = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vep.vcf"),
        v2t = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/reGenotype/vcf2table_MDS.py")
    conda: "./envs/vcf2csv.yml"
    threads: 4
    resources:
        mem_mb=20000,
        disk_mb=100000
    shell:
        """
        {BCFTOOLS} view -f 'PASS,clustered_events,multiallelic' {input.vcf} | python {input.v2t} simple --build GRCh37 -i /dev/stdin -o {output.csv}
        rm -rf /var/tmp/*
        """

rule annotate_mutect_cons:
    output:
        vcf = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vep.vcf")
    input:
        vcf = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vcf"),
        ref = S3.remote("{S3_BUCKET}/src/myelo_proliferation/mculen/ref/MIX_GRCh37-p13_GRCm38.p6-93.fa"),
        vep_data = S3.remote("{S3_BUCKET}/references/GRCh37-p13/test/VEP_95.tar.gz"),
    conda: "./envs/vep_new.yml"
    threads: 4
    resources:
        mem_mb=20000,
        disk_mb=100000
    shell:
        """
        tar -C {S3_BUCKET}/references/GRCh37-p13/test/ -xf {input.vep_data}
        {VEP} -i {input.vcf} --cache --cache_version 95 --dir_cache {S3_BUCKET}/references/GRCh37-p13/test/VEP_95 -fasta {input.ref} --merged --offline --vcf --everything -o {output.vcf}
        rm -rf /var/tmp/*
        """

rule mutect2_filter_normalize:
    output:
        S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.norm.vcf")
    input:
        S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.vcf")
        #ref = S3.remote("{S3_BUCKET}/references/GRCh37-p13/seq/GRCh37-p13.fa")
    conda: "./envs/bcftools.yml"
    threads: 4
    resources:
        mem_mb=10000,
        disk_mb=100000
    shell:
        """
        {BCFTOOLS} norm -m-both {input} > {output}
        rm -rf /var/tmp/*
        """

rule mutect2_filter:
     output:
         S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.filt.vcf")
     input:
         vcf = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.cons.vcf"),
     conda: "./envs/gatk4.yml"
     threads: 4
     resources:
         mem_mb=10000,
         disk_mb=100000
     shell:
         """
         gatk FilterMutectCalls -V {input.vcf} -O {output}
         rm -rf /var/tmp/*
         """

rule mutect2_cons:
    output:
        S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/variants/mutect2/{sample}.mutect2.cons.vcf")
    input:
        bam = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.second.bam"),
        ref = S3.remote("{S3_BUCKET}/src/myelo_proliferation/mculen/ref/MIX_GRCh37-p13_GRCm38.p6-93.fa"),
        fai = S3.remote("{S3_BUCKET}/src/myelo_proliferation/mculen/ref/MIX_GRCh37-p13_GRCm38.p6-93.fa.fai"),
        dct = S3.remote("{S3_BUCKET}/src/myelo_proliferation/mculen/ref/MIX_GRCh37-p13_GRCm38.p6-93.dict"),
        ivl = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.intervals"),
        bai = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.second.bam.bai")
    conda: "./envs/gatk4.yml"
    threads: 4
    resources:
         mem_mb=30000,
         disk_mb=100000
    shell:
        """
        gatk Mutect2 --reference {input.ref} --input {input.bam} --tumor-sample {wildcards.sample} --annotation StrandArtifact --min-base-quality-score 20 --intervals {input.ivl} --output {output}
        rm -rf /var/tmp/*
        """

##################
#### COVERAGE ####
##################

rule coverage_postprocess:
    input:
        txt = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/coverage/{sample}.PBcov.cons.txt"),
        rsc = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/scripts/coverage_stat_jana.R")
    output:
        S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/coverage/{sample}.perexon_stat.txt")
    conda: "./envs/erko.yml"
    threads: 4
    resources:
         mem_mb=10000,
         disk_mb=100000
    shell:
        """
        Rscript --vanilla {input.rsc} {input.txt}
        rm -rf /var/tmp/*
        """

rule coverage:
    input:
        bam = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.second.bam"),
        bed = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.bed")
    output:
        S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/coverage/{sample}.PBcov.cons.txt")
    conda: "./envs/bedtools.yml"
    threads: 4
    resources:
        mem_mb=80000,
        disk_mb=100000
    shell:
        """
        {BEDTOOLS} coverage -abam {input.bed} -b {input.bam} -d > {output}
        rm -rf /var/tmp/*
        """

#########################
#### QUALITY CONTROL ####
#########################

rule sedond_bam_QC:
    output:
        flagstat    = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/stats/second_bam_qc/{sample}.flagstat"),
        samstats    = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/stats/second_bam_qc/{sample}.samstats"),
        hs_metrics  = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/stats/second_bam_qc/{sample}.hs_metrics"),
        aln_metrics = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/stats/second_bam_qc/{sample}.aln_metrics")
    input:
        bam = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.second.bam"),
        jar = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/picard.jar"),
        ivl = S3.remote("{S3_BUCKET}/src/myelo_proliferation/archer/beds/jana_archer_unique_plus2nt.intervallist"),
        ref = S3.remote("{S3_BUCKET}/src/myelo_proliferation/mculen/ref/MIX_GRCh37-p13_GRCm38.p6-93.fa")
    conda: "./envs/samtools.yml"
    threads: 4
    resources:
        mem_mb=10000,
        disk_mb=100000
    shell:
        """
        {SAMTOOLS} flagstat {input.bam} > {output.flagstat}
        {SAMTOOLS} stats {input.bam} > {output.samstats}
        java -jar {input.jar} CollectHsMetrics I={input.bam} BAIT_INTERVALS={input.ivl} TARGET_INTERVALS={input.ivl} R={input.ref} O={output.hs_metrics}"
        java -jar {input.jar} CollectAlignmentSummaryMetrics I={input.bam} R={input.ref} O={output.hs_metrics}
        rm -rf /var/tmp/*
        """

########################
#### PREPROCESSING #####
########################

rule dedup:
    output:
        bam = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.second.bam"),
        bai = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.second.bam.bai")
    input:
        bam = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.first.bam"),
        bai = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.first.bam.bai")
    conda: "./envs/umi_tools.yml"
    log:
        stats = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.dedup")
    threads: 4
    resources:
        mem_mb=20000,
        disk_mb=100000
    shell:
           """
           {UMI_TOOLS} dedup -I {input} --paired --output-stats={log.stats} -S {output}
           {SAMTOOLS} index {output.bam}
           rm -rf /var/tmp/*
           """

rule align_bam:
    output:
        bam = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.first.bam"),
        bai = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.first.bam.bai")
    input:
        fwd = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R2.fastq"),
        ref_bwa = S3.remote(expand(S3_BUCKET + "/src/myelo_proliferation/mculen/ref/BWA/MIX_GRCh37-p13_GRCm38.p6-93{prip}", prip = ["", ".amb",".ann",".bwt",".pac",".sa"]))
    conda: "./envs/bwa.yml"
    log:
        run = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/bam/{sample}.log")
    threads: 4
    resources:
        mem_mb=60000,
        disk_mb=100000
    shell:
           """
           {BWA} mem -t 20 -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina' -v 1 {input.ref_bwa[0]} {input.fwd} {input.rev} 2>> {log.run} | {SAMTOOLS} view -Sb - | {SAMBAMBA} sort -o {output.bam} /dev/stdin 2>> {log.run}
           {SAMTOOLS} index {output.bam}
           rm -rf /var/tmp/*
           """

rule trim_adaptors:
    output:
        fwd = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.UMI.R2.fastq")
    input:
        fwd = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R2.fastq")
    conda: "./envs/cutadapt.yml"
    threads: 4
    resources:
        mem_mb=10000,
        disk_mb=100000
    shell: """
           {CUTADAPT} -g AACCGCCAGGAGT -m 50 -o {output.fwd} -p {output.rev} {input.fwd} {input.rev} > trim.out
           rm -rf /var/tmp/*
           """

rule UMI_extract:
    output:
        fwd = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R1.fastq"),
        rev = S3.remote("{S3_BUCKET}/data/mculen/{run}/samples/{sample}/fastq_trimmed/{sample}.trimmed.R2.fastq")
    input:
        unpack(getFasta)
    conda: "./envs/umi_tools.yml"
    threads: 4
    resources:
        mem_mb=10000,
        disk_mb=100000
    shell: """
           {UMI_TOOLS} extract -I {input.fwd} --bc-pattern=NNNNNNNN --read2-in={input.rev} --stdout={output.fwd} --read2-out={output.rev}
           rm -rf /var/tmp/*
           """
