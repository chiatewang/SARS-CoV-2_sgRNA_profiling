#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
## Title: Summarize intersecting ncsgRNAs
## Author: Chia-Te Wang
## Date: 2022-12
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
    #The purpose of this script is to process Venn Diagram intersecting csv file
    #Intersecting sgRNA and their relative expression will be write into tsv file
    #All_intersect_ncsg.tsv and All_intersect_Omicron_ncsg.tsv file will subsequently visualize using customized R script

# 1. Check arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--Variant_sets",required=True,
    help="e.g. -v 'BA.2, BA.5, BA.4, B.1.1.529, BA.1'")

args = parser.parse_args()
Variant_sets = args.Variant_sets

# 2. NcsgRNAs analysis (seperate each sample)

## import required packages
import os
import glob
import pandas as pd
cwd = os.getcwd()
intersect = pd.read_csv(cwd+"/"+"sgRNA_intersection_list.csv")
## In this filter part, users should key in "Variant sets" manually 
## (You can get the information from sgRNA_intersection_list.csv)
## e.g. 'BA.2, BA.5, BA.4, B.1.1.529, BA.1'
intersect_filter = intersect[intersect['Variant sets'] == Variant_sets]
intersect_filter = intersect_filter['Intersecting sgRNA'].values[0].split(', ')
path = glob.glob(cwd+"/raw_tsv/"+"*_sgRNA_nc_junction_summary.tsv")
#Ref. https://www.business-science.io/python/2021/09/21/python-read-csv.html

summary_col = []
summary_col = pd.DataFrame(summary_col, columns= ["lineage","j5","j3","count","startpos","productsize","JS","deletion size",
        "sample","RPM","chrom","type","name","start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'])
summary_col.to_csv("All_intersect_Omicron_ncsg.tsv",sep="\t",index=0)

## for loop parsing all tsv file in current directory
for i in path:
    RAW_TSV = pd.read_csv(i,sep="\t",names = ["j5","j3","count","startpos","productsize","JS","deletion size",
        "sample","RPM","chrom","type","name","start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'])
    prefix = i.split('/')
    variant_name = prefix[-1].split('_')
    RAW_TSV.insert(loc=0, column='lineage', value=variant_name[0])
    for k in intersect_filter:
        RAW_TSV1 = RAW_TSV[RAW_TSV['JS'] == k]
        RAW_TSV1.to_csv("All_intersect_Omicron_ncsg.tsv",sep="\t",mode="a+",index=0,header=0)

