# Tool Names??

## Tool Names??: SARS-CoV-2_sgRNA_profiling


### Description:

This repository and pipeline is for the work reported in the manuscript "Profiling subgenomic RNAs from different variants of SARS-CoV-2."

The purpose of this pipeline is to process Next Generataion Sequencing (NGS) data. 
This script maps the reads against the SARS-CoV-2 genome and generates junction coordinates from splice alignments; junction corrdinates next process using python script to classify sgRNAs and finally visualize using python script and R packages.

All python and bash scripts can be run independently. Execute any python script or bash script (without arguments) for detailed instructions on how to run them.


### Environments setup: 

1. This pipeline is running on Ubuntu 20.04.4 LTS (GNU/Linux 5.13.0-52-generic x86_64)
2. R (version 4.1)
3. python (version 3.6.10)


### Installation required tools:

1. Create conda environments and install some tools or packages
    $ conda create --name sgRNA python=3.6.10
    $ conda install -c bioconda samtools=1.9
    $ conda install -c anaconda pandas
    $ conda install -c conda-forge biopython
    $ conda install -c conda-forge matplotlib
2. Install Trimmomatic (version 0.39)   
    1. Download from http://www.usadellab.org/cms/?page=trimmomatic
    2. Add environment variable
        $ vi .bashrc
        add "export PATH={your install dir}/Trimmomatic-0.39" to .bashrc file
3. Install STAR mapping tool (version 2.7.9a)
    3.1 Download Source Code (tar.gz) from https://github.com/alexdobin/STAR/releases/tag/2.7.9a 
    3.2 $ tar -xzf 2.7.9a.tar.gz
        $ cd STAR-2.7.9a
    3.3 Compile
        $ cd STAR/source
        $ make STAR
    3.4 Add environment variable
        $ vi .bashrc
        add "export PATH={your install path}/STAR-2.7.9a/bin/Linux_x86_64" to .bashrc file
    3.5 Build genome index
        $ STAR \
        --runMode genomeGenerate -genomeDir ~/genomeDir_SARS2 \
        --runThreadN 4 \
        --genomeFastaFiles ~/SARS-CoV-2_sgRNA_profiling/resources/sars_cov2_NC_045512.2_genome.fasta \
        --sjdbGTFfeatureExon ~/SARS-CoV-2_sgRNA_profiling/resources/GCF_009858895.2_ASM985889v3_genomic.gff

### Usage:

1. Generate "sgRNA_junction_summary.tsv file"
    $ bash {your folder path}/sgRNA_pipeline.sh -d ~/GenomeDir_SARS2 -v {variant_name} -t {threads-numbers}
Note: sgRNA_pipeline.sh should be located at same directory with your fastq file
2. Data visualization
    $ R CMD BATCH venndiagram.r -i {"variant_name"_sgRNA_junction_summary.tsv}
