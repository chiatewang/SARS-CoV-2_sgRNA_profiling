## Subgenomic RNAs Profile ANalysis (SPAN): Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants.

### Description:
* A SgRNA Profile ANalysis (SPAN) is proposed for profiling various subgenomic RNAs (sgRNAs) and accompany with SPANviewer for data visualization. This pipeline is currently under submission in the manuscript titled **"Profiling of Noncanonical Subgenomic RNAs in SARS-CoV-2 Variants".**
* This pipeline first maps the reads to the SARS-CoV-2 reference genome and generates junction coordinates from splice alignments. The junction coordinates are then processed using a Python script to identify subgenomic RNAs and finally visualized using R Shiny. All scripts can be run independently.

### Requirements: 
1. Linux or MacOS
2. R (version 4.1)
3. python (version 3.8.10)

### Installation for SPAN:
Environments setup:

```
conda env create -f ./SPAN.yml -n SPAN
conda activate SPAN
```        

### Data visualization in our R shiny app (Plot_Profiling.R):
* Before using this app, please install R packages first.

```Rscript install.R```

Three pages we provided are as follows.

1. Venn Diagram page 
2. Junction Site Sashimi Plot page
3. Boxplot page

### Usage:  
1. Activate conda environment:
	```
	conda activate SPAN
	```
2. Please prepare fastq files in ~/SPAN/{variant_name} folder (e.g., use SRAToolkit fasterq-dump function to retrieve SRA reads)
3. Generate sample.txt:
	```
	python ./scripts/generate_sample_list.py
	```    
4. Execute this SPAN pipeline with snakemake:
	```
	snakemake -c1
	```
5. Data visualization:
	```
	Rscript ./R/SPANviewer.R
	```

### Authors @ Chang Gung University
* Yu-Nong Gong
* Chia-Te Wang