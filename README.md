# Subgenomic RNAs Profile ANalysis (SPAN)

## Subgenomic RNAs Profile ANalysis (SPAN): Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants.


### Description:

This repository and pipeline is for the work reported in the manuscript "Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants."

The purpose of this pipeline is to process Next Generation Sequencing (NGS) data. 
This script maps the reads against the SARS-CoV-2 genome and generates junction coordinates from splice alignments; junction coordinates next process using python script to identify subgenomic RNAs and finally visualize using R packages.

All scripts can be run independently. Execute any script for detailed instructions on how to run them.


### Hareware/software requirements: 

1. Linux and MacOS
2. R (version 4.1)
3. python (version 3.8.10)

### Installation:

1. Environments setup:

    $ cd {path_to_SPAN_folder}/SPAN
    $ conda env create -f ~/SPAN/SPAN.yml -n SPAN
    $ conda activate SPAN
    $ Rscript install.R
         
2. Build genome index:
    
    $ cd SPAN
    $ STAR --runMode genomeGenerate -genomeDir ./GenomeDir --runThreadN 4 --genomeFastaFiles ./resources/sars_cov2_NC_045512.2_genome.fasta --sjdbGTFfeatureExon ./resources GCF_009858895.2_ASM985889v3_genomic.gff
         
### Usage:  
1. Activate conda environment.
  
    $ conda activate SPAN
    
2. Prepare fastq files in ~/SPAN/{variant_name} folder
    
3. Generate "ncsgRNA_junction_summary.tsv"
  
    $ bash SPAN.sh -d GenomeDir -v {variant_name} -t {threads-numbers}  
     #### Example code:
    $ bash SPAN.sh -d GenomeDir -v B.1.1.529 -t 4  
    $ bash SPAN.sh -d GenomeDir -v BA.5 -t 4

4. Generate JSON file for Venn Diagram

    $ python output_to_json.py
    
5. Generate Venn Diagram and prepare "sgRNA_intersection_list.csv" from shiny app

    $ Rscript VennDiagram.R  
     #### Notes: VennDiagram visualizes using R shiny and users can control output using sidebar
    
6. Summarize intersecting ncsgRNAs from "sgRNA_intersection_list.csv" and "ncsgRNA_junction_summary.tsv"
    
    $ python Summarize_intersecting_ncsgRNAs.py -v {Variants_sets}  
     #### Notes: Variants_sets information can be found in "sgRNA_intersection_list.csv" and make sure to
     #### add " ' " at the begin and the end or you will get "error: unrecognized arguments:"
    
7. Data visualization  

    **VennDiagram.R** to generate venndiagram in the manuscript "Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants."  
    **SPAN_plot.R** to generate sashimi plot and ncsgRNAs expression plot in the manuscript "Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants."
    
    $ Rscript SPAN_plot.R -i sample_list.txt -f {lineage_name} -n {number_of_samples filter} -r {reads_per_million filter}  
     #### Notes: sample_list.txt which contains all lineage names except main lineage should be manually created by users.  
     #### Notes: -f choose the lineage in your comparsion that you want to be the main lineage,   
     #### -n filter number of samples greater than the number you type  
     #### (e.g. -n 4 means filter ncsgRNAs identified at least 5 different idependent samples)  
     #### Example Code:
    $ Rscript SPAN_plot.R -i sample_list.txt -f BA.5 -n 0 -r 10
