import pandas as pd
import os
import re
import argparse
import numpy

parser=argparse.ArgumentParser(description='Arguments to parse_bamreadcount.py')
parser.add_argument('--file', required=True)
parser.add_argument('--outname', required=True)


args=parser.parse_args()
filename=args.file
outname = args.outname

bamreadcount_files = list()

f = open(filename,'r')
for file in f.readlines():
    bamreadcount_files.append(file.rstrip())

# Create a template from the original BED file:
coverage = pd.read_csv('bp2bp_hg19.mtdna.bed',sep='\t')
# If necessary: coverage['ID'] = coverage['ID'].str.replace('MT_', 'chrM_')

reference_dict = dict()
sample_cols = list()

for i in bamreadcount_files:
    df = pd.read_csv(i,sep='\t', header= None)
    df['ID']  = df[[0, 1]].astype(str).apply(lambda x: '_'.join(x), axis=1)
    sample_name = i.split('.')[0].split('/')[-1]

    cov_dict = dict(zip(df['ID'], df[3]))
    A_dict = dict(zip(df['ID'], df[4]))
    C_dict = dict(zip(df['ID'], df[5]))
    G_dict = dict(zip(df['ID'], df[6]))
    T_dict = dict(zip(df['ID'], df[7]))
    reference_dict.update(dict(zip(df['ID'], df[2])))
    coverage[sample_name + '_cov'] = coverage['ID'].map(cov_dict)
    coverage[sample_name + '_A'] = coverage['ID'].map(A_dict) # This columns keep always the same order in a BmaReadCount output file. In any case, check it out!
    coverage[sample_name + '_C'] = coverage['ID'].map(C_dict)
    coverage[sample_name + '_G'] = coverage['ID'].map(G_dict)
    coverage[sample_name + '_T'] = coverage['ID'].map(T_dict)
    sample_cols.append([sample_name + '_cov',sample_name + '_A',sample_name + '_C',sample_name + '_G',sample_name + '_T'])#,sample_name' + _ACGT'


sampleid = list()

coverage['REFERENCE'] = coverage['ID'].map(reference_dict)

flat_sample_cols = [item for sublist in sample_cols for item in sublist]
cols_order = ['ID', 'CHROM','POS', 'REFERENCE','GENE'] + flat_sample_cols
coverage = coverage[cols_order]
coverage.to_csv(outname + '.BAMREADCOUNT_metrics.csv',sep='\t',index=None)
