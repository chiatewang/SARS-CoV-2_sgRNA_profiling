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
parser.add_argument("-n", "--FileNamePrefix",required=True)

args = parser.parse_args()
junction_path = args.junction_path
FileNamePrefix = args.FileNamePrefix

#-------------------------------------------------------------------------------------------------#

# 2. Get annotation file (You can change your annotation file here)
import pandas as pd
annotations = pd.read_csv('./SPAN/resources/SARS-CoV-2-annotations.gff', sep='\t',
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
## 3-1. Remove junction count <10 to avoid calculating error
junction_output = junction_output[junction_output['count'] >= 10 ]

#-------------------------------------------------------------------------------------------------#

## 4. Import reference (Wuhan strains)
import Bio
from Bio import SeqIO, Seq
covseq = str(next(SeqIO.parse('./SPAN/resources/sars_cov2_NC_045512.2_genome.fasta','fasta')).seq)

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

## 6. combine j5 and j3 into JS
junction_anno['JS'] = junction_anno['j5'].astype(str) + " - " + junction_anno['j3'].astype(str)
#junction_anno = junction_anno.groupby('JS').agg({'count': 'sum', 'j5': 'first', 'j3': 'first','startpos':'first','productsize':'first'})
# https://www.geeksforgeeks.org/ways-to-filter-pandas-dataframe-by-column-values/
junction_anno = junction_anno.sort_values(by=['count'], ascending=False)
## small deletion (potential deletion from template sequence)
junction_anno['deletion size'] = junction_anno['j3']-junction_anno['j5']
small_del = junction_anno[junction_anno['deletion size'] <= 1000 ]
## 7. remove small deletion
junction_anno = junction_anno[junction_anno['deletion size'] > 1000 ]

#-------------------------------------------------------------------------------------------------#

## 8. Classify canonical sgRNA
CANONICAL_START = 40 # 0-based, inc
CANONICAL_END = 80 # 0-based, noninc

is_cano = junction_anno['j5'].between(CANONICAL_START, CANONICAL_END-1)
ld_cano = junction_anno[is_cano]
ld_noncano = junction_anno[~is_cano]

can_prod = pd.merge(ld_cano, cdsanno, left_on='startpos', right_on='start')
can_prod = can_prod[can_prod['name'].notnull()]
can_prod['JS'] = can_prod['j5'].astype(str) + " - " + can_prod['j3'].astype(str)
#https://stackoverflow.com/questions/25333044/saving-dataframe-names-as-csv-file-names-in-pandas
can_prod.groupby('name').agg({'count': 'sum'}).sort_values(by='count', ascending=False).to_csv(FileNamePrefix + "canonical_sgRNA_count.tsv" , sep="\t")

#-------------------------------------------------------------------------------------------------#

## 9. Classify TRS-L noncanonical sgRNA
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

## 9-1. Test if startcodon is generated by the fusion
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
        #fusestarted = fusestarted.reindex(columns=['count','j5', 'j3' ,'startpos','productsize', 'deletion size', 'chorm' ,'type' ,'name' ,'start' ,'end','JS','evtype', 'evrelorf', 'evsize','startpos_byfusion']).head()
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
        #remaining = remaining.reindex(columns=['count','j5', 'j3' ,'startpos','productsize', 'deletion size', 'chorm' ,'type' ,'name' ,'start' ,'end','JS','evtype', 'evrelorf', 'evsize','startpos_byfusion']).head()
        #remaining.head()
else :
    print("no can_nonprod")


# 9-2. Classify TRS-L dependent noncanonical sgRNA (inframe or outframe)
if can_nonprod['count'].sum() > 0 :
    if remaining['count'].sum() == 0 :
        TRSL_inf=fusestarted.loc[fusestarted['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum()
        TRSL_outf=fusestarted.loc[fusestarted['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum()
    if fusestarted['count'].sum() == 0 :
        TRSL_inf=remaining.loc[remaining['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum()
        TRSL_outf=remaining.loc[remaining['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum()
    if remaining['count'].sum() != 0 and fusestarted['count'].sum() != 0 : 
        TRSL_inf=remaining.loc[remaining['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum() + fusestarted.loc[fusestarted['evtype'].isin(['N-term addition', 'N-term trunc']), 'count'].sum()
        TRSL_outf=remaining.loc[remaining['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum() + fusestarted.loc[fusestarted['evtype'].isin(['frameshift,mid', 'noORF','frameshift,5p']), 'count'].sum()
else :
    print("no can_nonprod")
    TRSL_inf = 0
    TRSL_outf = 0

#https://www.statology.org/pandas-sum-column-with-condition/

#-------------------------------------------------------------------------------------------------#

## 10. Classify TRS-L independent Non-canonical sgRNA
ld_noncano = pd.merge(ld_noncano, cdsanno, left_on='startpos', right_on='start', how='left')
ld_noncano = ld_noncano[ld_noncano['name'].isnull()]
ld_noncano['JS'] = ld_noncano['j5'].astype(str) + " - " + ld_noncano['j3'].astype(str)

delsize = ld_noncano['j3'] - ld_noncano['j5']
nc_distal = ld_noncano[delsize >= 5000].copy()
nc_proximal = ld_noncano[delsize < 5000].copy()
#nc_distal['count'].sum(), nc_proximal['count'].sum()

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

annotation = nc_distal.apply(annotate_orf_by_joint_position, axis=1)
annotation = annotation.apply(pd.Series)
annotation.columns = ['evtype', 'evrelorf', 'evsize']
nc_distal = pd.concat([
    nc_distal, annotation], axis=1)

nc_distal.groupby('evtype').agg({'count': 'sum'})

annotation = nc_proximal.apply(annotate_orf_by_joint_position, axis=1)
annotation = annotation.apply(pd.Series)
annotation.columns = ['evtype', 'evrelorf', 'evsize']
nc_proximal = pd.concat([
    nc_proximal, annotation], axis=1)

nc_proximal.groupby('evtype').agg({'count': 'sum'})

# 10-1. Classify TRS-L independent noncanonical sgRNA (inframe or outframe)

D_inf=nc_distal.loc[nc_distal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()
D_outf=nc_distal.loc[~nc_distal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()
P_inf=nc_proximal.loc[nc_proximal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()
P_outf=nc_proximal.loc[~nc_proximal['evtype'].isin(['in-frame-fusion-orf']), 'count'].sum()

#-------------------------------------------------------------------------------------------------#

# 11. Noncanonical sgRNA count summary

ncsg=['TRSL-inf','TRSL-outf','D-inf','D-outf','P-inf','P-outf']
ncsgcount=[TRSL_inf,TRSL_outf,D_inf,D_outf,P_inf,P_outf]

#https://www.delftstack.com/zh-tw/howto/python-pandas/pandas-create-dataframe-from-list/
ncsg_count = pd.DataFrame(list(zip(ncsg,ncsgcount)), columns = ['orf','count'])

ncsg_count.to_csv(FileNamePrefix + "ncsgRNA_count.tsv", sep="\t", mode='a+',index=False)

#-------------------------------------------------------------------------------------------------#

## 11. Generate sgRNA junction summary file
sgRNA_summary = []
#sgRNA_summary = pd.DataFrame(sgRNA_summary, columns= ['count','j5', 'j3' ,'startpos','productsize', 'deletion size', 'chorm' ,'type' ,'name' ,'start' ,'end','JS','evtype', 'evrelorf', 'evsize','startpos_byfusion'])
#sgRNA_summary.to_csv(FileNamePrefix + "sgRNA_junction_summary.tsv", sep="\t",index=False)
sgRNA_summary = pd.DataFrame(sgRNA_summary, columns= ['j5', 'j3' ,'count','startpos','productsize','JS' ,'deletion size', 'chorm' ,'type' ,'name' ,'start' ,'end','evtype', 'evrelorf', 'evsize','startpos_byfusion'])
sgRNA_summary.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t",index=False)
#https://blog.csdn.net/weixin_41529093/article/details/123411807
can_prod.to_csv(FileNamePrefix + "sgRNA_canonical_junction_summary.tsv", sep="\t",mode='a+',header=0,index=False)
if can_nonprod['count'].sum() > 0 :
    if fusestarted['count'].sum() == 0 :
        print("no fusestarted")
        remaining.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a+',header=0,index=False)
    else : 
        fusestarted.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t", mode='a+',header=0,index=False)
else :
    print("no can_nonprod")
nc_distal.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t",mode='a+',header=0,index=False)
nc_proximal.to_csv(FileNamePrefix + "sgRNA_nc_junction_summary.tsv", sep="\t",mode='a+',header=0,index=False)
small_del.to_csv(FileNamePrefix + "small_del.tsv", sep="\t",mode='a+',header=0,index=False)
noncanonical = pd.concat([ld_noncano,can_nonprod])
noncanonical['sgRNA ratio(%)'] = round(noncanonical["count"]/[junction_anno['count'].sum()]*100,4)
noncanonical = noncanonical.groupby('JS').agg({'count': 'sum', 'sgRNA ratio(%)':'sum','j5': 'first', 'j3': 'first','JS':'first'})
noncanonical = noncanonical.sort_values(by=['count'], ascending=False)                                  
noncanonical.head(50).to_csv(FileNamePrefix + "Top50.tsv", sep="\t",mode='a+',index=False)

#-------------------------------------------------------------------------------------------------#

## 12. Generate sashimi plot
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt
# define plot color
# color can be modify in "colorcode-rainbow.txt"
import pandas as pd
import matplotlib.colors as colors
orfcolors = pd.read_csv('/home/M1013001/ref/colorcode-rainbow.txt', sep=' ', names=['orf', 'color'], index_col=0)['color'].to_dict()
orfcolors

def common_tune(ax):
    ax.set_xlim(-.1, 29950)
    ax.set_ylim(-12, 12)

    for spname in 'top left right'.split():
        ax.spines[spname].set_visible(False)
    plt.setp(ax.get_yticklines(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)

ANCHOR_TOP = 1.5
ANCHOR_BOTTOM = -1.5
APEX_TOP = 9
APEX_BOTTOM = -9

def draw_bar(ax):
    for _, row in cdsanno.iterrows():
        if row['name'].startswith('nsp'):
            continue
        rect = Rectangle((row['start'], ANCHOR_BOTTOM), row['end'] - row['start'], ANCHOR_TOP-ANCHOR_BOTTOM, fc=orfcolors[row['name']])
        ax.add_patch(rect)

    rect = Rectangle((0, ANCHOR_BOTTOM), 265, ANCHOR_TOP-ANCHOR_BOTTOM, fc='#505050')
    ax.add_patch(rect)
    rect = Rectangle((265, ANCHOR_BOTTOM), 13441-265, ANCHOR_TOP-ANCHOR_BOTTOM, fc='#c5d9ea')
    ax.add_patch(rect)
    rect = Rectangle((13441, ANCHOR_BOTTOM), 21552-13441, ANCHOR_TOP-ANCHOR_BOTTOM, fc='#e8d2c4')
    ax.add_patch(rect)
    rect = Rectangle((29533, ANCHOR_BOTTOM), 29930-29533, ANCHOR_TOP-ANCHOR_BOTTOM, fc='#d0d0d0')
    ax.add_patch(rect)

def draw_junction(ax, start, end, top, height_coef=1, **kwds):
    distance = (end - start)

    height = max(1, height_coef * .015 * distance ** .6)
    curve = 3 * height ** .5
    pleft = start + distance / curve
    pmid = start + distance / 2
    pright = start + distance * (curve - 1) / curve
    anchor = ANCHOR_TOP if top else ANCHOR_BOTTOM
    apex = ANCHOR_TOP+height if top else ANCHOR_BOTTOM-height
    
    if 'edgecolor' in kwds:
        pass
    elif top:
        kwds['edgecolor'] = '#000000'
    else:
        kwds['edgecolor'] = '#DC143C'
    pp1 = mpatches.PathPatch(
        Path([(start, anchor), (pleft, apex), (pmid, apex), (pright, apex), (end, anchor)],
             [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.MOVETO]),
        fc="none", transform=ax.transData, **kwds)
    ax.add_patch(pp1)

import matplotlib.path as mpath
import matplotlib.patches as mpatches
Path = mpath.Path

ALPHA_MAP = [
    [1000, 1.0],
    [500, 0.5],
    [200, 0.3],
    [100, 0.2],
    [50, 0.1],
    [0, 0.01]
]
def calc_alpha(jnc):
    for mincnt, alpha in ALPHA_MAP:
        if jnc['count'] >= mincnt:
            return alpha
    return 0

fig, ax = plt.subplots(1, 1, figsize=(5.4, 2))

draw_bar(ax)

for _, jnc in nc_distal.iloc[:100].iterrows():
    alpha = calc_alpha(jnc)
    draw_junction(ax, jnc.j5, jnc.j3, top=('in-frame' not in jnc.evtype), alpha=alpha,
                  lw=.6)
    
common_tune(ax)
plt.savefig(FileNamePrefix +'plot-distal.pdf')

fig, ax = plt.subplots(1, 1, figsize=(5.4, 2))

draw_bar(ax)

for _, jnc in nc_proximal.iloc[:100].iterrows():
    alpha = calc_alpha(jnc)
    draw_junction(ax, jnc.j5, jnc.j3, top=('in-frame' not in jnc.evtype),
                  alpha=alpha, height_coef=5, lw=.6)
    
common_tune(ax)
plt.savefig(FileNamePrefix +'plot-proximal.pdf')

fig, ax = plt.subplots(1, 1, figsize=(5.4, 2))

draw_bar(ax)

for _, jnc in can_prod.iloc[:100].iterrows():
    alpha = calc_alpha(jnc)
    draw_junction(ax, jnc.j5, jnc.j3, top=True,
                  alpha=alpha, height_coef=1.2, lw=.6)
    
common_tune(ax)
plt.savefig(FileNamePrefix +'plot-canonical.pdf')

if can_nonprod['count'].sum() > 0 :
    if fusestarted['count'].sum() == 0 :
        rem_inframe = remaining[remaining['evtype'].isin(['N-term trunc'])]['j5 j3 count'.split()]
        rem_outframe = remaining[~remaining['evtype'].isin(['N-term trunc'])]['j5 j3 count'.split()]
        can_inframe = pd.concat([rem_inframe]).sort_values(by='count', ascending=False)
        can_outframe = pd.concat([rem_outframe]).sort_values(by='count', ascending=False)
        
        fig, ax = plt.subplots(1, 1, figsize=(5.4, 2))
        draw_bar(ax)
        
        for _, jnc in can_inframe.iloc[:100].iterrows():
            alpha = calc_alpha(jnc)
            draw_junction(ax, jnc.j5, jnc.j3, top=True,
                          alpha=alpha, height_coef=1, lw=.6)
            
        for _, jnc in can_outframe.iloc[:100].iterrows():
            alpha = calc_alpha(jnc)
            draw_junction(ax, jnc.j5, jnc.j3, top=False,
                          alpha=alpha, height_coef=1, lw=.6)
        common_tune(ax)
        plt.savefig(FileNamePrefix +'plot-TRSL.pdf')

    if remaining['count'].sum() == 0:
        fuse_inframe = fusestarted[fusestarted['evtype'].isin(['N-term addition'])]['j5 j3 count'.split()]
        fuse_outframe = fusestarted[~fusestarted['evtype'].isin(['N-term addition'])]['j5 j3 count'.split()]
        can_inframe = pd.concat([fuse_inframe]).sort_values(by='count', ascending=False)
        can_outframe = pd.concat([fuse_outframe]).sort_values(by='count', ascending=False)
        
        fig, ax = plt.subplots(1, 1, figsize=(5.4, 2))
        draw_bar(ax)
        
        for _, jnc in can_inframe.iloc[:100].iterrows():
            alpha = calc_alpha(jnc)
            draw_junction(ax, jnc.j5, jnc.j3, top=True,
                          alpha=alpha, height_coef=1, lw=.6)
            
        for _, jnc in can_outframe.iloc[:100].iterrows():
            alpha = calc_alpha(jnc)
            draw_junction(ax, jnc.j5, jnc.j3, top=False,
                          alpha=alpha, height_coef=1, lw=.6)
        common_tune(ax)
        plt.savefig(FileNamePrefix +'plot-TRSL.pdf')

    if remaining['count'].sum() != 0 and fusestarted['count'].sum() !=0 :
        fuse_inframe = fusestarted[fusestarted['evtype'].isin(['N-term addition'])]['j5 j3 count'.split()]
        fuse_outframe = fusestarted[~fusestarted['evtype'].isin(['N-term addition'])]['j5 j3 count'.split()]
        rem_inframe = remaining[remaining['evtype'].isin(['N-term trunc'])]['j5 j3 count'.split()]
        rem_outframe = remaining[~remaining['evtype'].isin(['N-term trunc'])]['j5 j3 count'.split()]
        can_inframe = pd.concat([fuse_inframe, rem_inframe]).sort_values(by='count', ascending=False)
        can_outframe = pd.concat([fuse_outframe, rem_outframe]).sort_values(by='count', ascending=False)
        
        fig, ax = plt.subplots(1, 1, figsize=(5.4, 2))
        draw_bar(ax)
        
        for _, jnc in can_inframe.iloc[:100].iterrows():
            alpha = calc_alpha(jnc)
            draw_junction(ax, jnc.j5, jnc.j3, top=True,
                          alpha=alpha, height_coef=1, lw=.6)
            
        for _, jnc in can_outframe.iloc[:100].iterrows():
            alpha = calc_alpha(jnc)
            draw_junction(ax, jnc.j5, jnc.j3, top=False,
                          alpha=alpha, height_coef=1, lw=.6)
        common_tune(ax)
        plt.savefig(FileNamePrefix +'plot-TRSL.pdf')   
        
else :
    print("no can_nonprod")

#-------------------------------------------------------------------------------------------------#
