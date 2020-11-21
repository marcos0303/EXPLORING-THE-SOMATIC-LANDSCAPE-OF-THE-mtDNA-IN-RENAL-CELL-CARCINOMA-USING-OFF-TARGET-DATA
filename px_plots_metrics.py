import plotly.express as px # pip install plotly https://plotly.com/python/plotly-express/
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler # pip install scikit-learn
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from scipy import stats

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
samples = list(dff['samples']) + list(dff['samples'])
mtcov = list(dff['tumor_mean_mtcov']) + list(dff['germline_mean_mtcov'])
totreads = list(dff['tumor_totalreads']) + list(dff['germline_totalreads'])
offtarget = list(dff['tumor_offtarget']) + list(dff['germline_offtarget'])
category = ['tumor']*len(dff) + ['germline']*len(dff)
gf = pd.DataFrame({'samples':samples,'Total reads':totreads,'mtDNA coverage':mtcov,'Sample type':category, 'Off target':offtarget})
gf['Total Off-target reads'] = gf['Total reads'] * gf['Off target']
tot_offtarget_reads = list(gf['Total Off-target reads'])
gf['Histological subtype'] = gf['samples'].map(dict(zip(list(db['sample']),list(db['histology_simplified']))))
gf['Primary/Mtx'] = gf['samples'].map(dict(zip(list(db['sample']),list(db['type_simplified']))))


for cat in list(set(list(gf['Sample type']))):
    cat_mtcov = list(gf['mtDNA coverage'].loc[gf['Sample type'] == cat])
    cat_offreads = list(gf['Total Off-target reads'].loc[gf['Sample type'] == cat])
    slope, intercept, r_value, p_value, std_err = stats.linregress(cat_offreads,cat_mtcov)
    print(cat)
    print('Slope:' + str(slope) + ', intercept:' + str(intercept)+ ', R:' + str(r_value)+ ', Pval:' + str(p_value)+ ', stderr:' + str(std_err))

#tumor
#Slope:1.3941953025203277e-05, intercept:-0.8269256090075139, R:0.5707113732284123, Pval:1.2155391401328011e-13, stderr:1.6953757057604007e-06
#germline
#Slope:2.0344470973947944e-06, intercept:-0.47704254845252514, R:0.767748873779222, Pval:7.538136797315453e-29, stderr:1.4350016295284426e-07


fig = px.scatter(gf, x='Total reads',y='mtDNA coverage', color='Sample type', marginal_x='box',marginal_y='violin', trendline="ols", template="simple_white", size='Off target', size_max=15, hover_name='samples') #, x='Total Off-target reads',log_x=True
fig.update_layout(
    legend=dict(
        x=0,
        y=1,
        xanchor='left',
        yanchor='top',
        traceorder="reversed",
        title_font_family="Courier",
        title="",
        font=dict(
            family="Courier",
            size=30,
            color="black"
        ),
        bgcolor="whitesmoke",
        bordercolor="Black",
        borderwidth=1
    )
)
fig.update_xaxes(tickfont=dict(size=25,color='black'),title_font=dict(size=35, family='Courier', color='black'))
fig.update_yaxes(tickfont=dict(size=25,color='black'),title_font=dict(size=35, family='Courier', color='black'))
fig.show()

## Same figure (OFF-line)


fig = px.scatter(gf, x='Total reads',y='mtDNA coverage', color='Sample type', marginal_x='box',marginal_y='violin', trendline="ols", template="simple_white", size='Off target', size_max=15, hover_name='samples') #, x='Total Off-target reads',log_x=True
fig.update_layout(
    legend=dict(
        x=0,
        y=1,
        xanchor='left',
        yanchor='top',
        traceorder="reversed",
        title_font_family="Courier",
        title="",
        font=dict(
            family="Courier",
            size=30,
            color="black"
        ),
        bgcolor="whitesmoke",
        bordercolor="Black",
        borderwidth=1
    )
)
fig.update_xaxes(tickfont=dict(size=25,color='black'),title_font=dict(size=35, family='Courier', color='black'))
fig.update_yaxes(tickfont=dict(size=25,color='black'),title_font=dict(size=35, family='Courier', color='black'))
plotly.offline.plot(fig, filename='/home/jlanillos/CNIO/IES Las Musas/try.png')


# Supp figure

fig = px.scatter(gf, x='Total reads',y='mtDNA coverage', color='Sample type', symbol='Histological subtype', marginal_x='box',marginal_y='box', template="simple_white", size='Off target', size_max=15, hover_name='samples') #, x='Total Off-target reads',log_x=True
fig.update_layout(
    legend=dict(
        x=0.8,
        y=1,
        xanchor='left',
        yanchor='top',
        traceorder="reversed",
        title_font_family="Courier",
        title="",
        font=dict(
            family="Courier",
            size=30,
            color="black"
        ),
        bgcolor="whitesmoke",
        bordercolor="Black",
        borderwidth=1
    )
)
fig.update_xaxes(tickfont=dict(size=25,color='black'),title_font=dict(size=35, family='Courier', color='black'))
fig.update_yaxes(tickfont=dict(size=25,color='black'),title_font=dict(size=35, family='Courier', color='black'))
fig.show()


#### Piechart

fig = px.sunburst(gf.loc[gf['Sample type'] == 'tumor'].groupby(['Histological subtype', 'Primary/Mtx']).count()['samples'].reset_index(), path=['Histological subtype', 'Primary/Mtx'], values='samples',names='Histological subtype',color='Histological subtype', hover_data=['Histological subtype','Primary/Mtx'])
colors = ['gold', 'mediumturquoise', 'darkorange', 'lightgreen']

fig.update_traces(textfont_size=20,marker=dict(colors=colors, line=dict(color='#000000', width=2)))
fig.update_layout(title={'text': "Histological subtypes and tumor location sites", 'x': 0.5, 'y': 0.97, 'xanchor' : 'center', 'yanchor': 'top', 'font' : dict(size=30, color="black")})

fig.show()
