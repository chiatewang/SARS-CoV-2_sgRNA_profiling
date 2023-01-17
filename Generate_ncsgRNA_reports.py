#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
## Title: Generate ncsgRNA reports
## Author: Chia-Te Wang
## Date: 2022-12
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
    #The purpose of this script is to process SPAN output tsv file
# 1. Check arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--Variant_sets",required=True,
    help="e.g. -v 'BA.2, BA.5, BA.4, B.1.1.529, BA.1'")

args = parser.parse_args()
Variant_sets = args.Variant_sets

## import required packages and working path
import os
import glob
import pandas as pd
cwd = os.getcwd()
#Ref. https://www.business-science.io/python/2021/09/21/python-read-csv.html

# 2. Noncanonical sgRNA summary
summary_col = []
summary_col = pd.DataFrame(summary_col, columns= ['lineage','JS','j5','j3','number_of_reads','RPM','number_of_samples','evtype'])
summary_col.to_csv("./supplemental_file_S2_Noncanonical_sgRNA_summary.csv",sep=",",index=0)
path = glob.glob(cwd+"/"+"raw_tsv/"+"*_sgRNA_nc_junction_summary.tsv")

## for loop parsing all tsv file in current directory
for i in path:
    RAW_TSV = pd.read_csv(i,sep="\t",names = ["j5","j3","count","startpos","productsize","JS","deletion size",
        "sample","RPM","chrom","type","name","start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'])
    TSV_summary = RAW_TSV.groupby('JS').agg({'j5':'first','j3':'first','count': 'sum','RPM':'sum','sample':'count','evtype':'first'}).sort_values(by='RPM', ascending=False).reset_index()
    TSV_summary.columns = ['JS','j5','j3','number_of_reads','RPM','number_of_samples','evtype']
    #TSV_summary = TSV_summary[TSV_summary['number_of_samples'] >= 5]
    prefix = i.split('/')
    variant_name = prefix[-1].split('_')
    TSV_summary.insert(loc=0, column='lineage', value=variant_name[0])
    TSV_summary.to_csv("supplemental_file_S2_Noncanonical_sgRNA_summary.csv",sep=",",mode="a+",index=0,header=0)

# 3. Print out canonical sgRNA count
summary_col = []
summary_col = pd.DataFrame(summary_col, columns= ['lineage','name', 'number_of_reads', 'RPM','number_of_samples'])
summary_col.to_csv("./supplemental_file_S3_Canonical_sgRNA_summary.csv",sep=",",index=0)
path = glob.glob(cwd+"/"+"raw_tsv/"+"*_canonical_sgRNA_count.tsv")

## for loop parsing all tsv file in current directory
for i in path:
    prefix = i.split('/')
    variant_name = prefix[-1].split('_')
    print(variant_name[0])
    CNO_TSV = pd.read_csv(i,sep="\t",names = ["name","count","RPM","sample"])
    CNO_summary = CNO_TSV.groupby('name').agg({'count': 'sum','RPM':'sum','sample':'count'}).sort_values(by='count', ascending=False).reset_index()
    CNO_summary.columns = ['name', 'number_of_reads', 'RPM','number_of_samples']
    CNO_summary.insert(loc=0, column='lineage', value=variant_name[0])
    CNO_summary.to_csv("supplemental_file_S3_Canonical_sgRNA_summary.csv",sep=",",mode="a+",index=0,header=0)

# 4. Omicron potential ncsgRNA analysis
intersect = pd.read_csv(cwd+"/"+"sgRNA_intersection_list.csv")
intersect_filter = intersect[intersect['Variant sets'] == Variant_sets]
intersect_filter = intersect_filter['Intersecting sgRNA'].values[0].split(', ')
path = glob.glob(cwd+"/"+"raw_tsv/"+"*_sgRNA_nc_junction_summary.tsv")

summary_col = []
summary_col = pd.DataFrame(summary_col, columns= ['lineage','JS','number_of_reads','RPM','number_of_samples','evtype'])
summary_col.to_csv("supplemental_file_S4_Omicron_intersect_ncsgRNA_summary.tsv",sep="\t",index=0)

## for loop parsing all tsv file in current directory
for i in path:
    RAW_TSV = pd.read_csv(i,sep="\t",names = ["j5","j3","count","startpos","productsize","JS","deletion size",
        "sample","RPM","chrom","type","name","start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'])
    TSV_summary = RAW_TSV.groupby('JS').agg({'count': 'sum','RPM':'sum','sample':'count','evtype':'first'}).sort_values(by='RPM', ascending=False).reset_index()
    TSV_summary.columns = ['JS','number_of_reads','RPM','number_of_samples','evtype']
    prefix = i.split('/')
    variant_name = prefix[-1].split('_')
    TSV_summary.insert(loc=0, column='lineage', value=variant_name[0])
    for k in intersect_filter:
        TSV_summary1 = TSV_summary[TSV_summary['JS'] == k]
        TSV_summary1.to_csv("supplemental_file_S4_Omicron_intersect_ncsgRNA_summary.tsv",sep="\t",mode="a+",index=0,header=0)
