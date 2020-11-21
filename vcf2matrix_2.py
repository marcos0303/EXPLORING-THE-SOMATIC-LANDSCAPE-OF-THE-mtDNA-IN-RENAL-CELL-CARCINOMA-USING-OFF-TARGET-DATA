import sys
import re
import argparse
import pandas as pd
import numpy as np
import numbers

parser=argparse.ArgumentParser(description='VCF file to transforme and output sample name')
parser.add_argument('--vcf', required=True)
parser.add_argument('--sample', required=True)
args=parser.parse_args()
filename=args.vcf
sample_name=args.sample
infile=open(filename,'r')
for l in infile:

#metadata
    if l.startswith('##INFO=<ID=CSQ,'): # get the info form the CSQ header
        CSQ_HEADER=re.compile('Format: (.*)">').search(l).group(1).split('|')
        continue
    if l.startswith('##'): continue     #other '##' lines unused
    if l.startswith('#CHROM'): #'#CHROM' is the header
        HEADER=l.strip('#').split()
        SAMPLES=HEADER[9:] # The sample name is on the 9th (and next) field of the VCF file header.
            # Create the Output HEADER (Filter: it's the one established after merging VCF files with bcftools merge (default); DP_allsamples (sum of the DP of all samples containing that variant (bftools merge -i DP:sum)
        OUTPUT_HEADER = ['ID','CHROM','POS','REF','ALT','FILTER'] + CSQ_HEADER #,'Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','cDNA_position','VARIANT_CLASS','HGNC_ID','Protein_position','Amino_acids','Codons','Existing_variation','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','CLIN_SIG','SIFT','PolyPhen','SOMATIC','PHENO'] # check the INFO header of your VCF for more ouput options
        SAMPLES_HEADER = []
        for sample in SAMPLES: SAMPLES_HEADER = SAMPLES_HEADER + ['_'.join([sample,'gt'])] + ['_'.join([sample,'ad'])] + ['_'.join([sample,'af'])]
        AF_samplenames = []
        for sample in SAMPLES: AF_samplenames = AF_samplenames + ['_'.join([sample,'af'])] # SAMPLE_NAMES_af to be used later during filtering
        OUTPUT= []
        continue

#variants

    # s will contain all fields for this line (including CSQ); Store all fields of a variant in a dictionary
    s=l.strip().split('\t')
    s=dict(zip(HEADER, s))

    # Add a new field containing each sample's name
    for k in SAMPLES: s[k]=s[k]

# Extract variant fields
    ## CSQ field
    csq_aux = [tuple(x.split('=')) for x in s['INFO'].split(';')]  # Ex: [('DP', '37'), ('ECNT', '1'), ('POP_AF', '5.000e-08'), ('P_CONTAM', '0.00'), ('P_GERMLINE', '-8.142e+00'), ('TLOD', '4.21'), ('CSQ', 'T|intron_v...]
    csq_aux = [y for y in csq_aux if len(y)>1]
    # SAMPLE_INFO_DP: the addition of all samples' sequencing depth
    SAMPLE_INFO_DP = int(dict(csq_aux)['DP'].split(',')[0])
    # SAMPLE_INFO_CSQ. useful to extract transcripts info. Ex: ['T|intron_variant|MODIFIER|STAG2|ENSG00000101972|Transcript|ENST00000218089|p...]
    SAMPLE_INFO_CSQ = dict(csq_aux)['CSQ'].split(',')
    # Parameters of interest to be included in the output file
    ID = '_'.join([s['CHROM'],s['POS'],s['REF'],s['ALT']]) # Variant ID containing the chrom, pos, ref and alt


# Extract information per sample
    SAMPLE_INFO_fields = []
    OUTPUT_general = [ID, s['CHROM'],s['POS'],s['REF'],s['ALT'],s['FILTER']]
    for SAMPLE in SAMPLES:
        SAMPLE_INFO = dict(zip(s['FORMAT'].split(':'),s[SAMPLE].split(':'))) # Ex: {'GT': '0/1', 'AD': '33,2', 'AF': '0.081', ... , 'SA_MAP_AF': '0.00,0.061,0.057', 'SA_POST_PROB': '0.011,0.019,0.970'}
        SAMPLE_INFO_fields = SAMPLE_INFO_fields + [SAMPLE_INFO['GT'], SAMPLE_INFO['AD'], SAMPLE_INFO['AF']] # SAMPLE_INFO['AD'][0], SAMPLE_INFO['AD'][1] are ref and alt allelic depths, respectively


# Preparing the output line and reading each annotated transcript for each variant. Sample info from CSQ added (see for loop below)
    for transcript in SAMPLE_INFO_CSQ: # There might be more than one transcript annotated, so they will be written in different lines
        CSQ = dict(zip(CSQ_HEADER,transcript.split('|')))
        CSQ= list(CSQ.values()) #[CSQ['Allele'],CSQ['Consequence'],CSQ['IMPACT'],CSQ['SYMBOL'],CSQ['Gene'],CSQ['Feature_type'],CSQ['Feature'],CSQ['BIOTYPE'],CSQ['EXON'],CSQ['cDNA_position'],CSQ['VARIANT_CLASS'],CSQ['HGNC_ID'],CSQ['Protein_position'],CSQ['Amino_acids'],CSQ['Codons'],CSQ['Existing_variation'],CSQ['gnomAD_AF'],CSQ['gnomAD_AFR_AF'],CSQ['gnomAD_AMR_AF'],CSQ['gnomAD_ASJ_AF'],CSQ['gnomAD_EAS_AF'],CSQ['gnomAD_FIN_AF'],CSQ['gnomAD_NFE_AF'],CSQ['gnomAD_OTH_AF'],CSQ['gnomAD_SAS_AF'],CSQ['CLIN_SIG'],CSQ['SIFT'],CSQ['PolyPhen'],CSQ['SOMATIC'],CSQ['PHENO']] # Create a list with the fiels from CSQ to be written in the output
        OUTPUT.append(OUTPUT_general + CSQ + SAMPLE_INFO_fields)


# OUTPUT list stored into pandas df
df = pd.DataFrame(OUTPUT,columns = OUTPUT_HEADER + SAMPLES_HEADER)

dff = df.loc[df['CANONICAL'] == 'YES'].copy()
dff.to_csv('MATRIX.canonical.csv',sep='\t',index = None)
df = dff.copy()
#############################################################################################################
########### ADD APPRIS, ALLELE FREQUENCIES ACROSS SAMPLES PER VARIANT, AND NUMBER OF SAMPLES PER VARIANT ####
#############################################################################################################
## Get _af columns with annotated AF, i.e., variants found
ll_af = [k for k in list(df.columns) if  k.endswith('_af')]
dff = df[ll_af].copy()
dff[dff[ll_af].applymap(lambda x: isinstance(x, numbers.Number)).all(1)]
dff = dff.fillna('.')
dff = dff.replace('.',0)
df['label'] = dff.apply(lambda x: ','.join(dff.columns[x.astype(bool)].tolist()), axis=1)
df['ALLELE_FREQUENCIES'] = dff[dff.columns[0:]].apply(lambda x: ','.join(x.dropna().astype(str)),axis=1)
df['ALLELE_FREQUENCIES'] = df['ALLELE_FREQUENCIES'].apply(lambda x: ','.join(list(filter(lambda y: y !='0', x.split(',')))))
df['ALLELE_FREQUENCIES'] = df['ALLELE_FREQUENCIES'].apply(lambda x: ','.join(list(filter(lambda y: y !='.', x.split(',')))))
df['ALLELE_FREQUENCIES'] = df['ALLELE_FREQUENCIES'].apply(lambda x: ','.join(list(filter(lambda y: y !='', x.split(',')))))


df['N_samples'] = np.array(list(map(len,df['label'].str.split(',').values)))
###########################################################################
ll_alt = [k for k in list(df.columns) if ('_ad' in k)]
#df_alt = df[ll_alt].copy()
df['AD_REF_ALT_allsamples'] = df[ll_alt].apply(lambda x: ';'.join(x.astype(str)),axis=1)

df['AD_REF_samples_vals'] = df['AD_REF_ALT_allsamples'].apply(lambda x: ','.join([i.split(',')[0] for i in x.split(';') if not '.' in i]))
df['AD_ALT_samples_vals'] = df['AD_REF_ALT_allsamples'].apply(lambda x: ','.join([i.split(',')[1] for i in x.split(';') if not '.' in i]))

df.to_csv('MATRIX.impact.canonical.csv',sep='\t',index = None)

################################################################################
########################### Create useful columns ##############################
#####################################  filter_choice() will save in a new column those samples passing HARD-FILTERING criteria (AF>=0.15  AND  #ALTreads > 2)
# It does not filter aou any variant (row) but just annotates new columns with those samples really passing our hard-filtering criteria
######################################################
def filter_choice(label, AF, ADREF, ADALT, minALT, minAF):
    AF_dict = dict(zip(label.split(','), AF.split(',')))
    ADREF_dict = dict(zip(label.split(','), ADREF.split(',')))
    ADALT_dict = dict(zip(label.split(','), ADALT.split(',')))
    AF_ADALT_dict = dict(zip(label.split(','),list(zip(AF.split(','), ADALT.split(',')))))

    AFfilter_dict = dict(filter(lambda elem: float(elem[1]) >= minAF,AF_dict.items()))
    ADALTfilter_dict = dict(filter(lambda elem: int(elem[1]) > minALT,ADALT_dict.items()))
    AF_ADALT_allfilters_dict = dict(filter(lambda elem: (float(elem[1][0]) >= minAF) and (int(elem[1][1]) > minALT),AF_ADALT_dict.items()))
    AF_allfilters_values = ','.join([value[0] for value in list(AF_ADALT_allfilters_dict.values())])
    ADALT_allfilters_values = ','.join([value[1] for value in list(AF_ADALT_allfilters_dict.values())])
    allfilters_samples = list(AF_ADALT_allfilters_dict.keys())
    ADREF_allfilters_values = ','.join([str(ADREF_dict[x]) for x in allfilters_samples])
    return ','.join(list(AFfilter_dict.keys())),','.join(list(AFfilter_dict.values())),','.join(list(ADALTfilter_dict.keys())),','.join(list(ADALTfilter_dict.values())),','.join(list(AF_ADALT_allfilters_dict.keys())),AF_allfilters_values,ADREF_allfilters_values,ADALT_allfilters_values


minALT = 1
minAF = 0.05
new_cols = ['AF_passed_samples', 'AF_passed_values', 'AD_ALT_passed_samples', 'AD_ALT_passed_values','All_filters_passed_samples', 'All_filters_passed_AFvalues','All_filters_passed_AD_REF_values','All_filters_passed_AD_ALT_values']
df[new_cols] = df.apply(lambda x: filter_choice(x['label'],x['ALLELE_FREQUENCIES'],x['AD_REF_samples_vals'],x['AD_ALT_samples_vals'], minALT, minAF), axis = 1, result_type="expand")

df['N_sample_passedFilters'] = df['All_filters_passed_samples'].apply(lambda x: len(x.split(',')))

new_cols = ['AF_passed_samples', 'AF_passed_values', 'AD_ALT_passed_samples', 'AD_ALT_passed_values','N_sample_passedFilters','All_filters_passed_samples', 'All_filters_passed_AFvalues','All_filters_passed_AD_REF_values','All_filters_passed_AD_ALT_values']

########################################################

# Prepare output matrix
OUTPUT_REHEADER = OUTPUT_HEADER + ['N_samples','label','ALLELE_FREQUENCIES','AD_REF_samples_vals','AD_ALT_samples_vals'] + new_cols +  SAMPLES_HEADER
df = df[OUTPUT_REHEADER].copy()
##################### SAVE RESULTS #####################
df = df.loc[df['All_filters_passed_samples'] != ''].copy()
df = df.loc[df['All_filters_passed_samples'] != 'NaN'].copy()
df.to_csv('matrix_' + sample_name + '.minALT' + str(minALT) + '.minAF' + str(minAF) + '.CANONICAL.csv',sep='\t',index = None)

