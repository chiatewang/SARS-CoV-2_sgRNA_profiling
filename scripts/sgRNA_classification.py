#!/usr/bin/env python
# coding: utf-8
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
## sgRNA_classification
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

# 1. Check arguments
import argparse
parser = argparse.ArgumentParser(description="""
            The purpose of this script is to classify junction-spanning reads from STAR 2.7.9a
            alignment result. 
            SJ_out.tab file provide information of JSRs that support sgRNAs
            Basically, this script echo 5' and 3' breakpoint and reads count of the discontinous 
            transcription 
            """)
parser.add_argument("-i", "--junction_path",required=True,
                    help="input STAR alignment SJ.out.tab")
parser.add_argument("-l", "--log",required=True,
                    help="input STAR alignment Log.final.out")
parser.add_argument("-s", "--sample",required=True)
parser.add_argument("-n", "--FileNamePrefix",required=True)

args = parser.parse_args()
junction_path = args.junction_path
log = args.log
sample = args.sample
FileNamePrefix = args.FileNamePrefix

#-------------------------------------------------------------------------------------------------#

# 2. Get annotation file (You can change your annotation file here)
import os
import glob
cwd = os.getcwd()
import pandas as pd
annotations = pd.read_csv(cwd+'/resources/SARS-CoV-2-annotations.gff', sep='\t',
                         names=['type', 'name', 'start', 'end'])
annotations['start'] = annotations['start']-1
#annotations['start'] -= 1
annotations = annotations.sort_values(by='start').reset_index(drop=True)

#-------------------------------------------------------------------------------------------------#

## 3. Import data (input STAR 2.7.9a splice alignment output Spliced Junction file)
junction_output = pd.read_csv(junction_path ,header=None, sep='\t')
junction_output = junction_output.iloc[:, [1,2,6]]
junction_output.columns = ['j5', 'j3','count']
## 3-1. Remove junction count <10 or count >0 to avoid calculating error
junction_output = junction_output[junction_output['count'] > 0 ]
#junction_output = junction_output[junction_output['count'] >= 10 ]

#-------------------------------------------------------------------------------------------------#

## 4. Import reference (Wuhan strains)
from Bio import SeqIO

for seq in SeqIO.parse(cwd+"/"+'/resources/sars_cov2_NC_045512.2_genome.fasta','fasta'):
    seq_id = seq.id
    refseq = str(seq.seq)

#-------------------------------------------------------------------------------------------------#

## 5. Combine j5 and j3 into JS
junction_output['JS'] = junction_output['j5'].astype(str) + " - " + junction_output['j3'].astype(str)
junction_output = junction_output.sort_values(by=['count'], ascending=False)
# https://www.geeksforgeeks.org/ways-to-filter-pandas-dataframe-by-column-values/

#-------------------------------------------------------------------------------------------------#

## 6. Remove small deletion (potential deletion from template sequence)
junction_output['deletion size'] = junction_output['j3']-junction_output['j5']
small_del = junction_output[junction_output['deletion size'] <= 1000 ] 
small_del.to_csv(FileNamePrefix + "small_del.tsv", sep="\t",mode='a',header=0,index=False)
junction_output = junction_output[junction_output['deletion size'] > 1000 ]
junction_output = junction_output.assign(sample=sample)

#-------------------------------------------------------------------------------------------------#

## 7. Normalized with mapped reads per million
Log = pd.read_csv(log,header=None, sep='\t')
Mapped_reads = float(Log.iat[7,1])
#http://justimchung.blogspot.com/2018/06/pandas-dataframe.html
junction_output['sgRPM'] = junction_output['count'].astype(float)*1000000/Mapped_reads
junction_output = junction_output.reset_index(drop=True)

#-------------------------------------------------------------------------------------------------#

## 8. Find "ATG" in reference sequence as start codon
startposition = []
for i in junction_output.index:
    seq_5 = refseq[:(junction_output['j5'][i])]
    seq_3 = refseq[(junction_output['j3'][i]):]
    seq_fuse = seq_5 + seq_3
    
    if 'ATG' in seq_fuse:
        startpos_fuse = seq_fuse.find('ATG')
        if startpos_fuse < (junction_output['j5'][i]):
            startpos = startpos_fuse
        else:
            startpos = startpos_fuse - (junction_output['j5'][i]) + (junction_output['j3'][i])
    else:
        startpos = 0
    startposition.append(startpos)
startposition = {'startpos':startposition}
startposition=pd.DataFrame(startposition)

junction_anno=pd.concat([junction_output,startposition],axis=1)

#-------------------------------------------------------------------------------------------------#

## 9. Classify canonical sgRNA
trs_L_start = 40 
trs_L_end = 80

trs_L_d = junction_anno['j5'].between(trs_L_start, trs_L_end)
sg_trs_L = junction_anno[trs_L_d]
sg_trs_L_ind = junction_anno[~trs_L_d]
can_sg = pd.merge(sg_trs_L, annotations, left_on='startpos', right_on='start')
can_sg = can_sg[can_sg['name'].notnull()]

can_sg.groupby('name').agg({'count': 'sum','sgRPM':'sum','sample':'first'}).sort_values(by='count', ascending=False).to_csv(FileNamePrefix + "canonical_sgRNA_count.tsv",sep="\t",header=0,mode='a')
can_sg.to_csv(FileNamePrefix + "sgRNA_canonical_junction_summary.tsv", sep="\t",mode='a',header=0,index=False)

#-------------------------------------------------------------------------------------------------#

## 10. Classify noncanonical sgRNA
trs_L_ncsg = pd.merge(sg_trs_L, annotations, left_on='startpos', right_on='start', how='left')
trs_L_ncsg = trs_L_ncsg[trs_L_ncsg['name'].isnull()]
if trs_L_ncsg['count'].sum() > 0 :
    trs_L_ncsg = trs_L_ncsg.iloc[:,:8]
    trs_L_ncsg.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False) 
else : 
    print("no trs_L_ncsg")
if sg_trs_L_ind['count'].sum() > 0 :
    sg_trs_L_ind.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False) 
else :
    print("no sg_trs_L_ind")

#-------------------------------------------------------------------------------------------------#
