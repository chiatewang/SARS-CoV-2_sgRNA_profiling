#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
## output_to_json
## Author: Chia-Te Wang
## Date:2022-12
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
    # The purpose of this script is to process Subgenomic RNA ANalysis output tsv file using python pandas and json.
    # Noncanonicial subgenomic RNA junction site will be summarized and converted into JSON format
    # JSON file will subsequently visualize using customized R script

## 1. import packages and set variables
import os
import glob
import pandas as pd
cwd = os.getcwd()
path = glob.glob(cwd+"/"+"raw_tsv/"+"*_sgRNA_nc_junction_summary.tsv")
#Ref. https://www.business-science.io/python/2021/09/21/python-read-csv.html
dict_lineage = {}

# Main Function
## 2. import tsv and summarize JS with "number of reads" and "number of samples"
## for loop to parsing all ncsgRNA tsv file at current directory
for i in path:
    RAW_TSV = pd.read_csv(i,sep="\t",names = ["j5","j3","count","startpos","productsize","JS","deletion size",
        "sample","RPM","chrom","type","name","start","end",'evtype', 'evrelorf', 'evsize','startpos_byfusion'])
    TSV_summary = RAW_TSV.groupby('JS').agg({'count': 'sum','RPM':'sum','sample':'count',"j5":"first",'evtype':'first'}).sort_values(by='sample', ascending=False).reset_index()
    TSV_summary.columns = ['JS', 'number_of_reads', 'RPM','number_of_samples','j5','evtype']
    ## Filter datasets 
    # Users can determine the filter options
    # Example script (filter datasets that only echo ncsgRNAs number of samples >= 5)
    #TSV_summary = TSV_summary[TSV_summary['number_of_samples'] > 4 ]

    ## write dataframe to dictionary
    #Ref : https://stackoverflow.com/questions/26716616/convert-a-pandas-dataframe-to-a-dictionary
    dict_info = TSV_summary.to_dict('list')

    prefix = i.split('/')
    variant_name = prefix[-1].split('_')
    print(variant_name[0])
    #https://steam.oxxostudio.tw/category/python/basic/dictionary.html
    #https://ithelp.ithome.com.tw/articles/10220160
    dict_lineage[variant_name[0]]= dict_info

## 3. compile different lineages into one json file
# Generate json for venn diagram input 
import json
with open("sgRNA_list.json", "w") as outfile:
    json.dump(dict_lineage, outfile,indent=len(path))