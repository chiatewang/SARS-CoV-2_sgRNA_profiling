#!/bin/bash

## Defining usage and setting input

usage() {
        echo "
        1. Reference Mapping
         Mapping short reads to a reference genome using Spliced Transcripts 
         Alignment to a Reference (STAR) aligner (Dobin et al.)
         Parameters are adapted from Kim et al's settings.
        2. Junction Site Identification
         Spliced Junction output must be parse by sgRNA_classification.py adapted from Kim et al. to 
         get subgenomic RNA (sgRNA) classification and product prediction
         Filtering 5' end >=40 and <=80 (instead of 55~85 mentioned in Kim et al.) to define 
         TRS-L dependent sgRNA. 
         Classification of sgRNA (Canonical, Noncanonical)
         Normolization by calculating total sgRNA counts
        
        USAGE: $0

        Required:
        -d <INDEX>
            Path to the STAR index directory. 
        -v <VARIANTS>
            Target variants that you want to profile
        Optional:
        -t <THREADS> [1]
            Number of computing threads available.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] 
then
        echo "Input Error"
        echo "Usage: bash SPAN.sh -d GenomeDir -v {variant_name} -t {threads-numbers}"
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
mkdir raw_tsv
path="$PWD"
cd $VARIANTS
date +"%y-%m-%d-%T">> "$VARIANTS"_time.log
START_total=`date +"%s"`

# Data preprocessing (Trim Illumina Universal Adapters and low quality reads)
# Capture Acc number and set a variable
for i in `ls *_1.fastq`
    do
        Forward=${i##*/}
        Forward=${Forward%.*}
        Reverse=${Forward%_*}

        echo "1. Data preprocessing."
        #cut illumina universal adapters
        trimmomatic PE -phred33 \
        "$Forward".fastq "$Reverse"_2.fastq \
        "$Forward"_paired.fastq "$Forward"_unpaired.fastq \
        "$Reverse"_2_paired.fastq "$Reverse"_2_unpaired.fastq \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
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
        --genomeDir "$path"/GenomeDir/ \
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
        --alignMatesGapMax 1000000 \
        --limitBAMsortRAM 110000000000
        
        rm *pair*
        cd "$path"
        echo "3. Junction Site Identification."
        python "$path"/sgRNA_classification.py \
        -i `echo "$path"/"$VARIANTS"/"$Reverse"_SJ.out.tab` \
        -l `echo "$path"/"$VARIANTS"/"$Reverse"_Log.final.out` \
        -s "$Reverse" \
        -n "$VARIANTS"_
done

echo "Script Finished"
date +"%y-%m-%d-%T">> "$VARIANTS"_time.log
END_total=`date +"%s"`
let RUNNING_total=($END_total-$START_total)
echo "script running time = " $RUNNING_total " seconds" >> "$VARIANTS"_time.log
mv "$VARIANTS"* "$path"/raw_tsv/

