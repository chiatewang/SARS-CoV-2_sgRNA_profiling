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
                         names=['chrom', 'type', 'name', 'start', 'end', '_x', '_y', '_z'])
annotations['start'] -= 1

# cds with Leader sequence
cdsanno = annotations[(annotations['type'] == 'CDS') &
                      (~annotations['name'].isin(['ORF10', 'nsp1', 'frameshift'])) &
                      (~annotations['name'].apply(lambda x: x.startswith('nsp')))]
cdsanno = pd.concat([
    cdsanno,
    pd.DataFrame([
        pd.Series(['chrSCV', 'CDS', 'ORF1a', 265, 13467, 0.0, '.', '.'], index=cdsanno.columns),
        pd.Series(['chrSCV', 'CDS', 'ORF1b', 13467, 21552, 0.0, '.', '.'], index=cdsanno.columns),
    ])
]).sort_values(by='start').reset_index(drop=True).iloc[:, :5]

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
import Bio
from Bio import SeqIO, Seq
covseq = str(next(SeqIO.parse(cwd+"/"+'/resources/sars_cov2_NC_045512.2_genome.fasta','fasta')).seq)

#-------------------------------------------------------------------------------------------------#

## 5. Find "ATG" in reference sequence as start codon
junction_anno = []
for _, row in junction_output.iterrows():
    seqleft = covseq[:row.j5]
    seqright = covseq[row.j3:]
    seqrecomb = seqleft + seqright

    if 'ATG' in seqrecomb:
        # calculate start codon positions
        startpos_recomb = seqrecomb.find('ATG')
        if startpos_recomb < row.j5:
            startpos = startpos_recomb
        else:
            startpos = startpos_recomb - row.j5 + row.j3

        # calculate protein product size
        orfcandseq = seqrecomb[startpos_recomb:]
        if len(orfcandseq) % 3 > 0:
            orfcandseq += 'N' * (3 - len(orfcandseq) % 3)

        translation = Seq.Seq(orfcandseq).translate().split('*')[0]
        #if row['name'] == 'j2':
        #    print(Seq(orfcandseq).translate())
        product_size = len(translation)
    else:
        startpos_recomb = -1
        startpos = -1
        product_size = 0
        translation = ''

    junction_anno.append([row['j5'], row['j3'], row['count'],
                          startpos, product_size])
junction_anno = pd.DataFrame(junction_anno, columns=['j5', 'j3', 'count', 'startpos', 'productsize'])

#-------------------------------------------------------------------------------------------------#

## 6. Combine j5 and j3 into JS
junction_anno['JS'] = junction_anno['j5'].astype(str) + " - " + junction_anno['j3'].astype(str)
junction_anno = junction_anno.sort_values(by=['count'], ascending=False)
# https://www.geeksforgeeks.org/ways-to-filter-pandas-dataframe-by-column-values/

## 7. Remove small deletion (potential deletion from template sequence)
junction_anno['deletion size'] = junction_anno['j3']-junction_anno['j5']
small_del = junction_anno[junction_anno['deletion size'] <= 1000 ] 
junction_anno = junction_anno[junction_anno['deletion size'] > 1000 ]
junction_anno = junction_anno.assign(sample=sample)
#-------------------------------------------------------------------------------------------------#

## 8. Normalized with mapped reads per million
Log = pd.read_csv(log,header=None, sep='\t')
Mapped_reads = float(Log.iat[7,1])
#http://justimchung.blogspot.com/2018/06/pandas-dataframe.html
junction_anno['sgRPM'] = junction_anno['count'].astype(float)*1000000/Mapped_reads

#-------------------------------------------------------------------------------------------------#

## 9. Classify canonical sgRNA
CANONICAL_START = 40 
CANONICAL_END = 80 

is_cano = junction_anno['j5'].between(CANONICAL_START, CANONICAL_END-1)
ld_cano = junction_anno[is_cano]
ld_noncano = junction_anno[~is_cano]

can_prod = pd.merge(ld_cano, cdsanno, left_on='startpos', right_on='start')
can_prod = can_prod[can_prod['name'].notnull()]
can_prod['JS'] = can_prod['j5'].astype(str) + " - " + can_prod['j3'].astype(str)
#https://stackoverflow.com/questions/25333044/saving-dataframe-names-as-csv-file-names-in-pandas
can_prod.groupby('name').agg({'count': 'sum','sgRPM':'sum','sample':'first'}).sort_values(\
    by='count', ascending=False).to_csv(FileNamePrefix + "canonical_sgRNA_count.tsv" , \
    sep="\t",header=0,mode='a')
can_prod.to_csv(FileNamePrefix + "sgRNA_canonical_junction_summary.tsv", sep="\t",mode='a',header=0,index=False)

#-------------------------------------------------------------------------------------------------#

## 10. Classify TRS-L noncanonical sgRNA
can_nonprod = pd.merge(ld_cano, cdsanno, left_on='startpos', right_on='start', how='left')
can_nonprod = can_nonprod[can_nonprod['name'].isnull()]
can_nonprod['JS'] = can_nonprod['j5'].astype(str) + " - " + can_nonprod['j3'].astype(str)

def measure_translation_length(startpos):
    seq = covseq[startpos:]
    if len(seq) % 3 > 0:
        seq += 'N' * (3 - (len(seq) % 3))
    return len(Seq.Seq(seq).translate().split('*')[0])

def annotate_orf(startpos):
    ovl_cds = (
        cdsanno[(cdsanno['start'] <= startpos) &
                (startpos < cdsanno['end'])])

    if len(ovl_cds) >= 2:
        ovl_cds = ovl_cds[~ovl_cds['name'].apply(lambda x: x.startswith('nsp') or x == 'ORF7b')]
        if len(ovl_cds) != 1:
            raise ValueError("Unresolvable conflict: " + ' '.join(ovl_cds['name']))
    elif len(ovl_cds) < 1:
        # Downstream
        nextcdses = cdsanno[cdsanno['start'] >= startpos]
        if len(nextcdses) == 0:
            return 'noORF', None, None
        nextcds = nextcdses.iloc[0]
        offset = nextcds['start'] - startpos
        if offset % 3 == 0:
            return 'N-term addition', nextcds['name'], offset // 3
        else:
            return 'frameshift,5p', nextcds['name'], offset
    
    # Single-hit within a known CDS
    ovl_cds = ovl_cds.iloc[0]
    offset = startpos - ovl_cds['start']
    if offset % 3 == 0:
        return 'N-term trunc', ovl_cds['name'], offset // 3
    
    neworflength = measure_translation_length(startpos)
    return 'frameshift,mid', ovl_cds['name'], neworflength

    # noORF
    # N-term addition
    # N-term trunc
    # frameshift,mid
    # frameshift,5p


## 10-1. Test if startcodon is generated by the fusion
def test_startcodon_fusiogenesis(row):
    fuse = covseq[int(row['j5'])-2:int(row['j5'])] + covseq[int(row['j3']):int(row['j3'])+2]
    if 'ATG' in fuse:
        position = fuse.index('ATG') # 0 or 1
        return 2-position
    else:
        return -1

if can_nonprod['count'].sum() == 0 :
    print("no can_nonprod")
else :
    fusestart_pos = can_nonprod.apply(test_startcodon_fusiogenesis, axis=1)
    fusestart_pos
    fusestarted = can_nonprod[fusestart_pos >= 0].copy()
    fusestarted['startpos_byfusion'] = fusestarted['j3'] - fusestart_pos[fusestart_pos >= 0]
    remaining = can_nonprod[fusestart_pos < 0]

if can_nonprod['count'].sum() > 0 :
    if fusestarted['count'].sum() == 0 :
        print("no fusestarted")
    else : 
        fusestarted_annotation = fusestarted['startpos_byfusion'].apply(annotate_orf)
        fusestarted_annotation = fusestarted_annotation.apply(pd.Series)
        fusestarted_annotation.columns = ['evtype', 'evrelorf', 'evsize']
        fusestarted = pd.concat([fusestarted,fusestarted_annotation], axis=1)
        #https://www.codeleading.com/article/24433444392/
else :
    print("no can_nonprod")

if can_nonprod['count'].sum() > 0 :
    if remaining['count'].sum() == 0 :
        print("no remaining")
    else : 
        remaining_annotation = remaining['startpos'].apply(annotate_orf)
        remaining_annotation = remaining_annotation.apply(pd.Series)
        remaining_annotation.columns = ['evtype', 'evrelorf', 'evsize']
        remaining = pd.concat([remaining,remaining_annotation], axis=1)
else :
    print("no can_nonprod")


# 10-2. Classify TRS-L dependent noncanonical sgRNA (inframe or outframe)
if can_nonprod['count'].sum() > 0 :
    if remaining['count'].sum() == 0 :
        TRSL_inf=fusestarted.loc[fusestarted['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum()
        TRSL_outf=fusestarted.loc[fusestarted['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum()
        fusestarted.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False)
    if fusestarted['count'].sum() == 0 :
        TRSL_inf=remaining.loc[remaining['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum()
        TRSL_outf=remaining.loc[remaining['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum()
        remaining.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False)
    if remaining['count'].sum() != 0 and fusestarted['count'].sum() != 0 : 
        TRSL_inf=remaining.loc[remaining['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum() + fusestarted.loc[fusestarted['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum()
        TRSL_outf=remaining.loc[remaining['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum() + fusestarted.loc[fusestarted['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum()
        remaining.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False)
        fusestarted.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False)
else :
    print("no can_nonprod")
    TRSL_inf = 0
    TRSL_outf = 0

#https://www.statology.org/pandas-sum-column-with-condition/

#-------------------------------------------------------------------------------------------------#

## 11. Classify TRS-L independent Non-canonical sgRNA
ld_noncano = pd.merge(ld_noncano, cdsanno, left_on='startpos', right_on='start', how='left')
ld_noncano = ld_noncano[ld_noncano['name'].isnull()]
ld_noncano['JS'] = ld_noncano['j5'].astype(str) + " - " + ld_noncano['j3'].astype(str)

delsize = ld_noncano['j3'] - ld_noncano['j5']
nc_distal = ld_noncano[delsize >= 5000].copy()
nc_proximal = ld_noncano[delsize < 5000].copy()

def measure_translation_length_by_seq(seq):
    if len(seq) % 3 > 0:
        seq += 'N' * (3 - (len(seq) % 3))
    return len(Seq.Seq(seq).translate().split('*')[0])

def annotate_orf_by_joint_position(row):
    j5, j3 = row['j5'], row['j3']
    j5_1 = j5 - 1 # last nucleotide in the 5' fragment

    ovl_cds5 = (
        cdsanno[(cdsanno['start'] <= j5_1) & (j5_1 < cdsanno['end'])])

    if len(ovl_cds5) >= 2:
        ovl_cds5 = ovl_cds5[~ovl_cds5['name'].apply(lambda x: x.startswith('nsp') or x == 'ORF7b')]
        if len(ovl_cds5) != 1:
            raise ValueError("Unresolvable conflict: {} ".format(j5_1) + ' '.join(ovl_cds5['name']))
    elif len(ovl_cds5) < 1:
        # j5 is in the middle of non-coding region. Use the downstream of j3
        nextcdses = cdsanno[cdsanno['start'] >= j3]
        if len(nextcdses) == 0:
            return 'no-orf', None, None
        nextcds = nextcdses.iloc[0]
        return 'may-translate', nextcds['name'], nextcds['start'] - j3
    
    # Single-hit within a known CDS
    ovl_cds5 = ovl_cds5.iloc[0]
    offset5 = j5 - ovl_cds5['start']
    seqleft = covseq[ovl_cds5['start']:j5]
    seqright = covseq[j3:]
    frame5 = len(seqleft) % 3

    # Find 3' position identity
    ovl_cds3 = (
        cdsanno[(cdsanno['start'] <= j3) & (j3 < cdsanno['end'])])
    if len(ovl_cds3) >= 2:
        ovl_cds3 = ovl_cds3[~ovl_cds3['name'].apply(lambda x: x.startswith('nsp') or x == 'ORF7b')]
        if len(ovl_cds3) != 1:
            raise ValueError("Unresolvable conflict: {} ".format(j3) + ' '.join(ovl_cds3['name']))
    elif len(ovl_cds3) < 1:
        # j5 is in ORF, j3 is not in ORF.
        product_len = measure_translation_length_by_seq(seqleft + seqright)

        nextcdses = cdsanno[cdsanno['start'] >= j3]
        if len(nextcdses) > 0:
            nextcds = nextcdses.iloc[0]
            nextstartoffset = nextcds['start'] - j3
            if product_len > (offset5 + nextstartoffset) / 3:
                return 'mid-utr-inserted-fusion', '{}-{}'.format(ovl_cds5['name'], nextcds['name']), (offset5/3, nextstartoffset/3)

        return 'C-term-trunc-and-extension-by-UTR', ovl_cds5['name'], (offset5/3, product_len)

    # Single-hit within a known CDS for 3' site too.
    ovl_cds3 = ovl_cds3.iloc[0]
    offset3 = j3 - ovl_cds3['start']
    remaining3 = ovl_cds3['end'] - j3
    frame3 = offset3 % 3
    
    if frame5 == frame3:
        return 'in-frame-fusion-orf', '{}-{}'.format(ovl_cds5['name'], ovl_cds3['name']), (offset5/3, remaining3/3)

    product_len = measure_translation_length_by_seq(seqleft + seqright)
    return 'frameshift-fusion-orf', '{}-{}'.format(ovl_cds5['name'], ovl_cds3['name']), (offset5/3, product_len)

annotate_orf_by_joint_position(nc_distal.iloc[1])

#no-orf
#may-translate
#mid-utr-inserted-fusion
#C-term-trunc-and-extension-by-UTR
#in-frame-fusion-orf
#frameshift-fusion-orf

## 11. Classify TRS-L independent noncanonical sgRNA (inframe or outframe)
if ld_noncano['count'].sum() > 0 :
    if nc_distal['count'].sum() == 0 :
        print("no nc_distal")
        D_inf = 0
        D_outf = 0
    else : 
        annotation = nc_distal.apply(annotate_orf_by_joint_position, axis=1)
        annotation = annotation.apply(pd.Series)
        annotation.columns = ['evtype', 'evrelorf', 'evsize']
        nc_distal = pd.concat([nc_distal, annotation], axis=1)
        D_inf=nc_distal.loc[nc_distal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()
        D_outf=nc_distal.loc[~nc_distal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()
        nc_distal.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False)      
else :
    print("no ld_noncano")

if ld_noncano['count'].sum() > 0 :
    if nc_proximal['count'].sum() == 0 :
        print("no nc_proximal")
        P_inf = 0
        P_outf = 0
    else : 
        annotation = nc_proximal.apply(annotate_orf_by_joint_position, axis=1)
        annotation = annotation.apply(pd.Series)
        annotation.columns = ['evtype', 'evrelorf', 'evsize']
        nc_proximal = pd.concat([nc_proximal, annotation], axis=1)
        P_inf=nc_proximal.loc[nc_proximal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()
        P_outf=nc_proximal.loc[~nc_proximal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()
        nc_proximal.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a',header=0,index=False)
else :
    print("no ld_noncano")

#-------------------------------------------------------------------------------------------------#

# 12. Noncanonical sgRNA count summary

ncsg=['TRSL-inf','TRSL-outf','D-inf','D-outf','P-inf','P-outf']
ncsgcount=[TRSL_inf,TRSL_outf,D_inf,D_outf,P_inf,P_outf]

#https://www.delftstack.com/zh-tw/howto/python-pandas/pandas-create-dataframe-from-list/
ncsg_count = pd.DataFrame(list(zip(ncsg,ncsgcount)), columns = ['orf','count'])
ncsg_count = ncsg_count.assign(sample=sample)
ncsg_count.to_csv(FileNamePrefix + "ncsgRNA_count.tsv",sep="\t",header=0,mode='a',index=False)

#-------------------------------------------------------------------------------------------------#

small_del.to_csv(FileNamePrefix + "small_del.tsv", sep="\t",mode='a',header=0,index=False)

#-------------------------------------------------------------------------------------------------#
