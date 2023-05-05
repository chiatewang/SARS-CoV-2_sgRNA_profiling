##########################################################
# 2022/12/09
# Chia-Te Wang
# Subgenomic RNA ANalysis Pipeline (SPAN) with Snakemake
##########################################################
# Snakemake frequently asked questions
# https://blog.csdn.net/yangl7/article/details/119249760

import pandas as pd
# ===== read in the sample list csv and convert into a dictionary structure
'''
{'Acc_number': {
  'fq1': '/Path/to/fastq/variant_name/fastq1',
  'fq2': '/Path/to/fastq/variant_name/fastq2',
  'variant_name': 'Omicron'}, ...
}
'''

def sample_list_csv_to_dict(sample_csv_path):
    df = pd.read_csv(sample_csv_path, sep='\t')
    df = df.dropna()
    df = df.set_index('id').T.to_dict('dict')
    return df

sample_list = sample_list_csv_to_dict('sample.txt')

def get_input_fastq_path(sample_wildcard, r):
   # a function getting fastq paths for specific sample id, r shall be either 1 or 2 for paired-end
    fastq_path = sample_list[sample_wildcard]['fq'+str(r)]
    return fastq_path

def assign_variant(sample_wildcard):
    variant = sample_list[sample_wildcard]['variant_name']
    return variant

# ========

rule all:
    input:
        expand("fastq/{variant}/{sample}_{R}.fastq", sample=list(sample_list.keys()), R=["1","2"], variant=glob_wildcards("{*}")),
        expand("results/cleanfastq/{sample}_1_paired.fastq", sample=list(sample_list.keys())),
        expand("results/cleanfastq/{sample}_1_unpaired.fastq", sample=list(sample_list.keys())),
        expand("results/cleanfastq/{sample}_2_paired.fastq", sample=list(sample_list.keys())),
        expand("results/cleanfastq/{sample}_2_unpaired.fastq", sample=list(sample_list.keys())),
        expand("{sample}_Aligned.sortedByCoord.out.bam", sample=list(sample_list.keys())),
        expand("{sample}_SJ.out.tab", sample=list(sample_list.keys())),
        expand("{sample}_Log.final.out", sample=list(sample_list.keys())),
        expand("{sample}_sgRNA_nc_junction_summary.tsv", sample=list(sample_list.keys())),
        expand("{sample}_sgRNA_canonical_junction_summary.tsv", sample=list(sample_list.keys())),
        "sgRNA_list.json"

rule data_preprocessing:
    input:
        r1 = lambda wildcards: get_input_fastq_path(wildcards.sample, 1),
        r2 = lambda wildcards: get_input_fastq_path(wildcards.sample, 2)
    output:
        o1 = "results/cleanfastq/{sample}_1_paired.fastq",
        o2 = "results/cleanfastq/{sample}_1_unpaired.fastq",
        o3 = "results/cleanfastq/{sample}_2_paired.fastq",
        o4 = "results/cleanfastq/{sample}_2_unpaired.fastq",
        l1 = "results/cleanfastq/{sample}_trim.log"
    shell:
        "trimmomatic PE -phred33 {input.r1} {input.r2} \
        {output.o1} {output.o2} {output.o3} {output.o4} \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 \
        SLIDINGWINDOW:4:15 \
        TRAILING:3 \
        MINLEN:36 \
       -trimlog {output.l1}"

rule generate_index:
    input:
        ref = "resources/sars_cov2_NC_045512.2_genome.fasta",
        gtf = "resources/GCF_009858895.2_ASM985889v3_genomic.gff"
    output:
        directory("GenomeDir_SARS2")
    shell:
        "STAR --runMode genomeGenerate \
        -genomeDir {output} \
        --runThreadN 4 \
        --genomeFastaFiles {input.ref} \
        --sjdbGTFfeatureExon {input.gtf} &&"
        "mv GenomeDir GenomeDir_SARS2"

rule map_reads:
    input:
        ref = "GenomeDir_SARS2",
        r1 = "results/cleanfastq/{sample}_1_paired.fastq",
        r2 = "results/cleanfastq/{sample}_2_paired.fastq"
    output:
        "{sample}_Aligned.sortedByCoord.out.bam",
        "{sample}_SJ.out.tab",
        "{sample}_Log.final.out"
    shell:
        "STAR \
        --runMode alignReads \
        --runThreadN 4 \
        --genomeDir {input.ref} \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix {wildcards.sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --outSJfilterOverhangMin 12 12 12 12 \
        --outSJfilterCountUniqueMin 1 1 1 1 \
        --outSJfilterCountTotalMin 1 1 1 1 \
        --outSJfilterDistToOtherSJmin 0 0 0 0 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --scoreGapNoncan -4 \
        --scoreGapATAC -4 \
        --chimOutType WithinBAM HardClip \
        --chimScoreJunctionNonGTAG 0 \
        --alignSJstitchMismatchNmax -1 -1 -1 -1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --limitBAMsortRAM 110000000000 "

rule junction_site_identification:
    input:
        i1 = "{sample}_SJ.out.tab",
        i2 = "{sample}_Log.final.out"
    params:
        prefix = lambda wildcards: assign_variant(wildcards.sample)
    output:
        o1 = "{sample}_sgRNA_nc_junction_summary.tsv",
        o2 = "{sample}_sgRNA_canonical_junction_summary.tsv"
    shell:
        "python scripts/sgRNA_classification.py -i {input.i1} \
        -l {input.i2} \
        -s {wildcards.sample}\
        -n {wildcards.sample}_ &&"
        "cat {output.o1} >> {params.prefix}_merge_sgRNA_nc_junction_summary.tsv &&"
        "cat {output.o2} >> {params.prefix}_merge_sgRNA_canonical_junction_summary.tsv "

rule output_to_json:
    output:
        "sgRNA_list.json"
    script:
        "scripts/output_to_json.py"