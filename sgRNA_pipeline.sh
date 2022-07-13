#!/bin/bash

## Defining usage and setting input

usage() {
        echo "
        Notice : fastq file must at same directory with this shell script
	The purpose of this shell script : 
        1. Reference Mapping
         Mapping short reads to a reference genome using Spliced Transcripts 
         Alignment to a Reference (STAR) aligner (Dobin et al.)
         Parameters are adapted from Kim et al's settings and is desired 
         to limit penalties for non-canonical read junctions.
        2. Junction Site Identification
         Spliced Junction output must be parse by sgRNA_classification.py adapted from Kim et al. to 
         get sibgenomic RNA (sgRNA) prediction
         Filtering 5' end >=40 and <=80 (instead of 55~85 mentioned in Kim et al.) to define 
         TRS-L dependent sgRNA. 
         Classification of sgRNA (Canonical, TRS-L dependent Novel, non-canonical distal 
         non-canonical proximal)
         Calculating sgRNA counts
        3. sgRNA Percentage Plot
         Generating boxplot to visualize each types of sgRNA
	
        USAGE: $0

        Required:
        -i <INPUT_FASTQ>
            Type fastq to parsing all .fastq file in current folder.
        -d <INDEX>
            Path to the star index directory. 
        Optional:
        -t <THREADS> [1]
            Number of computing threads available.
        "
}

#If less than 2 options are input, show usage and exit script.
if [ $# -le 2 ] 
then
        echo "Input Error"
        echo "Usage: bash sgRNA_pipeline.sh -d genomeDir -v variants (-t threads numbers)"
        exit 1
fi

#Setting input
while getopts d:v:t: option ; do
        case "${option}"
        in
                d) INDEX=${OPTARG};;
                v) VARIANTS=${OPTARG};;
                t) THREADS=${OPTARG};;
        esac
done

# Set defaults
THREADS=${THREADS:-1}

# Run cd-hit
echo "
INPUTS -

INDEX: $INDEX
VARIANTS: $VARIANTS
THREADS: $THREADS
"

echo "Starting script."
date +"%y-%m-%d-%T">> "$VARIANTS"_time.log
START_total=`date +"%s"`

# Data preprocessing (Trim Illumina Universal Adapters)
# Capture Acc number and set a variable
for i in `ls *_1.fastq`
do
        Forward=${i##*/}
        Forward=${Forward%.*}
        Reverse=${Forward%_*}

        echo "1. Data preprocessing."
        #cut illumina universal adapters
        java -jar /home/M1013001/lib/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
        "$Forward".fastq "$Reverse"_2.fastq \
        "$Forward"_paired.fastq "$Forward"_unpaired.fastq \
        "$Reverse"_2_paired.fastq "$Reverse"_2_unpaired.fastq \
        ILLUMINACLIP:/home/M1013001/lib/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
        #https://ithelp.ithome.com.tw/articles/10219913
        #https://www.biostars.org/p/478810/

        echo "2. Reference Mapping."
        
        STAR \
        --runMode alignReads \
        --runThreadN $THREADS \
        --genomeDir $INDEX \
        --readFilesIn "$Forward"_paired.fastq "$Reverse"_2_paired.fastq \
        --outFileNamePrefix $Reverse"_" \
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
        --alignMatesGapMax 1000000
        #--limitBAMsortRAM 8488872000
        # https://github.com/alexdobin/STAR/issues/870
        # For not enough memory for BAM sorting by using --outBAMsortingBinsN default to 200 or even more

        #"too short" literally means "alignment too short". This could either happen for normal-length 
        #read which are not mapping well for a read that was short (over-trimmed) before mapping,
        #The average trimmed length before mapping is reported is "average input length".
        #https://github.com/alexdobin/STAR/issues/577
        #https://github.com/alexdobin/STAR/issues/164
        #fasterq-dump *.sra to process paired-end reads

        # Generate index and coverage file
        samtools index `echo ./*bam`
        samtools depth -a `echo ./*bam` > $Reverse"_coverage".tsv

        cat ./"$Reverse"_SJ.out.tab >> "$VARIANTS"_SJ.out.tab
done

        # GenomeDir
        #STAR \
        #--runMode genomeGenerate -genomeDir /home/M1013001/genomeDir_SARS2 \
        #--runThreadN 4 \
        #--genomeFastaFiles /home/M1013001/ref/sars_cov2_NC_045512.2_genome.fasta \
        #--sjdbGTFfeatureExon /home/M1013001/GCF_009858895.2_ASM985889v3_genomic.gff

    ## STAR hyperparameter tuning
        # Apply Grid search concept ?
        # First to minimize file size
        # remove human genome
        # remove whole virus genome
        # remove canonical sgRNA
        # Second to pick out novel sgRNA that may alter by mapping parameter

echo "3. Junction Site Identification."

python /home/M1013001/workspace/script/sgRNA_classification.py \
-i `echo ./"$VARIANTS"_SJ.out.tab` \
-n "$VARIANTS"_

echo "4. sgRNA Percentage Plot"

python /home/M1013001/workspace/script/sgRNA_percentage.py \
-c `echo ./"$VARIANTS"_canonical_sgRNA_count.tsv` \
-s `echo ./"$VARIANTS"_nc_sgRNA_count.tsv` \
-n "$VARIANTS"_

echo "Script Finished"
date +"%y-%m-%d-%T">> "$VARIANTS"_time.log
END_total=`date +"%s"`
let RUNNING_total=($END_total-$START_total)
echo "script running time = " $RUNNING_total " seconds" >> "$VARIANTS"_time.log
