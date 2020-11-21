import plotly.express as px # pip install plotly https://plotly.com/python/plotly-express/
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler # pip install scikit-learn
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from scipy import stats
import matplotlib.gridspec as gridspec


mtgenes = pd.read_csv('/home/jlanillos/CNIO/IES Las Musas/mtDNA_genes.csv',sep='\t') # from https://genome.ucsc.edu/cgi-bin/hgTables

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


tumor['mean_mtdna_cov'] = tumor[[i for i in list(tumor.columns) if '_cov' in i]].mean(axis = 1)
#tumor['mean_mtdna_cov'].plot.line(color = 'cornflowerblue', linewidth = 2)
#germ['mean_mtdna_cov'] = germ[[i for i in list(germ.columns) if '_cov' in i]];lines = germ['mean_mtdna_cov'].plot.line(color = 'darkorange', linewidth = 2)
germ['mean_mtdna_cov'] = germ[[i for i in list(germ.columns) if '_cov' in i]].mean(axis = 1);
#lines = germ['mean_mtdna_cov'].plot.line(color = 'darkorange', linewidth = 2)


r = pd.DataFrame()
r['germline_mean_mtdna_cov'] =  germ['mean_mtdna_cov']
r['tumor_mean_mtdna_cov'] =  tumor['mean_mtdna_cov']
t = pd.DataFrame()
t['Tumor samples'] = r['tumor_mean_mtdna_cov']
t['Blood samples'] = r['germline_mean_mtdna_cov']
rr = t.copy()

#gf.groupby(['Histological subtype'])['mtDNA coverage'].mean()



fig = plt.figure(figsize=(10, 26))
gs = gridspec.GridSpec(2, 1, height_ratios=[1,3])
gs.update(wspace=0.085, hspace=0.01) # set the spacing between axes
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])
rr['Blood samples'].plot.line(color = 'darkorange', linewidth = 2, ax = ax2)
rr['Tumor samples'].plot.line(color = 'cornflowerblue', linewidth = 2, ax = ax2)
mtgenes[['DIFF']].T.plot.barh(width = 0.05, stacked=True, legend=None, ax = ax1)
mtgenes['textcoor'] = mtgenes['txStart'] + mtgenes['DIFF']/2
COORTEXT = list(mtgenes['textcoor'].loc[mtgenes['DIFF'] > 300])
TEXT = list(mtgenes['name2'].loc[mtgenes['DIFF'] > 300])
for i,j in list(zip(TEXT,COORTEXT)):
    ax1.text(j,0.05,i, rotation=45)

ax1.set_title('mtDNA genes (>300bp)')
ax1.axis('off')
ax2.set_xlabel('mtDNA genomic position')
ax2.set_ylabel('Mean mtDNA coverage')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.show()
