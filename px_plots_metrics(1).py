import plotly.express as px # pip install plotly https://plotly.com/python/plotly-express/
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler # pip install scikit-learn
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe

#germline_file_bamreadcount='/mnt/64716603-5b56-4f9a-b195-c11560647a3a/Projects/mtDNA_RCC_project_ieslasmusas/data/BAMREADCOUNT/germline.BAMREADCOUNT_metrics.onlyCOV.csv'
#tumor_file_bamreadcount='/mnt/64716603-5b56-4f9a-b195-c11560647a3a/Projects/mtDNA_RCC_project_ieslasmusas/data/BAMREADCOUNT/tumors.BAMREADCOUNT_metrics.onlyCOV.csv'
#rccdb='/mnt/64716603-5b56-4f9a-b195-c11560647a3a/Projects/mtDNA_RCC_project_ieslasmusas/data/RCC_samples_mtDNA_db_PAIRED_available_20201019.csv'
#germline_hsmetrics='/mnt/64716603-5b56-4f9a-b195-c11560647a3a/Projects/mtDNA_RCC_project_ieslasmusas/hsmetrics_germline/germline.HSALL.csv'
#tumor_hsmetrics='/mnt/64716603-5b56-4f9a-b195-c11560647a3a/Projects/mtDNA_RCC_project_ieslasmusas/hsmetrics_tumor/tumor.HSALL.csv'

germline_file_bamreadcount='/home/jlanillos/CNIO/IES Las Musas/data/germline.BAMREADCOUNT_metrics.onlyCOV.csv'
tumor_file_bamreadcount='/home/jlanillos/CNIO/IES Las Musas/data/tumors.BAMREADCOUNT_metrics.onlyCOV.csv'
germ = pd.read_csv(germline_file_bamreadcount,sep='\t');tumor = pd.read_csv(tumor_file_bamreadcount,sep='\t')
samples = list(tumor[[i for i in list(tumor.columns) if '_cov' in i]].mean(axis=0).index); values = list(tumor[[i for i in list(tumor.columns) if '_cov' in i]].mean(axis=0).values)
d = {'samples': samples, 'tumor_mean_mtcov':values}
dff = pd.DataFrame(d)
samples = list(germ[[i for i in list(germ.columns) if '_cov' in i]].mean(axis=0).index); values = list(germ[[i for i in list(germ.columns) if '_cov' in i]].mean(axis=0).values); d = dict(zip(samples,values))
dff['germline_mean_mtcov'] = dff['samples'].map(d)
dff['samples'] = dff['samples'].apply(lambda x: x.rstrip('_cov'))


germline_hsmetrics='/home/jlanillos/CNIO/IES Las Musas/data/germline.HSALL.csv'
tumor_hsmetrics='/home/jlanillos/CNIO/IES Las Musas/data/tumor.HSALL.csv'
hsallgerm = pd.read_csv(germline_hsmetrics,sep='\t')
hsalltumor = pd.read_csv(tumor_hsmetrics,sep='\t')
hsalltumor['tumorID'] = hsalltumor['BAIT_SET'].apply(lambda x: x.split('_')[0])
hsallgerm['germlineID'] = hsallgerm['BAIT_SET'].apply(lambda x: x.split('_')[0])


rccdb='/home/jlanillos/CNIO/IES Las Musas/data/RCC_samples_mtDNA_db_PAIRED_available_20201019.csv'
db = pd.read_csv(rccdb,sep='\t')
dff['tumorID'] = dff['samples'].map(dict(zip(list(db['sample']),list(db['tumor']))))
dff['germlineID'] = dff['samples'].map(dict(zip(list(db['sample']),list(db['blood']))))
dff['Cov20x'] = [1 if x >= 20 else 0 for x in dff['tumor_mean_mtcov']]
dff['tumor_totalreads'] = dff['tumorID'].map(dict(zip(list(hsalltumor['tumorID']),list(hsalltumor['TOTAL_READS']))))
dff['germline_totalreads'] = dff['germlineID'].map(dict(zip(list(hsallgerm['germlineID']),list(hsallgerm['TOTAL_READS']))))
dff['tumor_offtarget'] = dff['tumorID'].map(dict(zip(list(hsalltumor['tumorID']),list(hsalltumor['PCT_OFF_BAIT']))))
dff['germline_offtarget'] = dff['germlineID'].map(dict(zip(list(hsallgerm['germlineID']),list(hsallgerm['PCT_OFF_BAIT']))))
dff['tumor_bait_coverage'] = dff['tumorID'].map(dict(zip(list(hsalltumor['tumorID']),list(hsalltumor['MEAN_BAIT_COVERAGE']))))
dff['germline_bait_coverage'] = dff['germlineID'].map(dict(zip(list(hsallgerm['germlineID']),list(hsallgerm['MEAN_BAIT_COVERAGE']))))
dff['tumor_mtdnareads'] = dff['tumor_mean_mtcov']*16579/100
dff['germline_mtdnareads'] = dff['germline_mean_mtcov']*16579/100
dff['tumor_mtdnareads_totalreads_pct'] = dff['tumor_mtdnareads'] / dff['tumor_totalreads']
dff['germline_mtdnareads_totalreads_pct'] = dff['germline_mtdnareads'] / dff['germline_totalreads']
dff['norm_tumor_totalreads'] = dff['tumor_totalreads'] / np.max(dff['tumor_totalreads'])
dff['norm_germline_totalreads'] = dff['germline_totalreads'] / np.max(dff['germline_totalreads'])
dff['norm_tumor_mtdnareads'] = dff['tumor_mtdnareads'] / np.max(dff['tumor_mtdnareads'])
dff['norm_germline_mtdnareads'] = dff['germline_mtdnareads'] / np.max(dff['germline_mtdnareads'])
dff['pctnorm_tumor_readsratio'] = dff['norm_tumor_mtdnareads'] / dff['norm_tumor_totalreads']
dff['pctnorm_germline_readsratio'] = dff['norm_germline_mtdnareads'] / dff['norm_germline_totalreads']



## Plotly express plots
mtcov = list(dff['tumor_mean_mtcov']) + list(dff['germline_mean_mtcov'])
totreads = list(dff['tumor_totalreads']) + list(dff['germline_totalreads'])
offtarget = list(dff['tumor_offtarget']) + list(dff['germline_offtarget'])
category = ['tumor']*len(dff) + ['germline']*len(dff)
gf = pd.DataFrame({'totreads':totreads,'mtcov':mtcov,'category':category, 'offtarget':offtarget})
gf['Total Off-target reads'] = gf['totreads'] * gf['offtarget']



fig = px.scatter(gf, x='totreads',y='mtcov', color='category', marginal_x='histogram',marginal_y='histogram')
fig.show()
