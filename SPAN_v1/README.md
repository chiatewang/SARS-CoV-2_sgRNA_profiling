# Subgenomic RNAs Profile ANalysis (SPAN)

## Subgenomic RNAs Profile ANalysis (SPAN): Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants.


### Description:

This repository and pipeline is for the work reported in the manuscript "Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants."

The purpose of this pipeline is to process Next Generation Sequencing (NGS) data. 
This script maps the reads against the SARS-CoV-2 genome and generates junction coordinates from splice alignments; junction coordinates next process using python script to identify subgenomic RNAs and finally visualize using R packages.

All scripts can be run independently. Execute any script for detailed instructions on how to run them.


### Hareware/software requirements: 

1. Linux or MacOS
2. R (version 4.1)
3. python (version 3.8.10)

### Installation:

1. Environments setup:

    $ cd {path_to_SPAN_folder}/SPAN_v1
    $ conda env create -f ~/SPAN_v1/SPAN.yml -n SPAN
    $ conda activate SPAN
    $ Rscript install.R
         
2. Build genome index:
    
    $ cd SPAN_v1
    $ STAR --runMode genomeGenerate -genomeDir ./GenomeDir --runThreadN 4 --genomeFastaFiles ./resources/sars_cov2_NC_045512.2_genome.fasta --sjdbGTFfeatureExon ./resources GCF_009858895.2_ASM985889v3_genomic.gff
         
### Usage:  
1. Activate conda environment.
  
    $ conda activate SPAN
    
2. Prepare fastq files in ~/SPAN/{variant_name} folder
    
3. Generate "ncsgRNA_junction_summary.tsv"
  
    $ bash SPAN.sh -d GenomeDir -v {variant_name} -t {threads-numbers}  
     #### Example code:
    $ bash SPAN.sh -d GenomeDir -v BA.1 -t 4  
    $ bash SPAN.sh -d GenomeDir -v BA.2 -t 4

4. Generate JSON file for Data visualization

    $ python output_to_json.py
    
5. Data visualization using R shiny

    $ Rscript Plot_Profiling.R  
    