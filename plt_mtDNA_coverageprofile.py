import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
vectorsnames = ['vectorgermline','vectortumor']
files = ['germline.BAMREADCOUNT_metrics.csv',
'tumors.BAMREADCOUNT_metrics.csv']

filename_dict = dict(zip(vectorsnames,files))
vector_dict = dict()
for i in filename_dict.keys():
    filename = filename_dict[i]
    df = pd.read_csv(filename,sep='\t')
    cov_cols = [x for x in list(df.columns) if '_cov' in x]
    dff = pd.DataFrame()
    dff[cov_cols] = df[cov_cols].copy()
    normdff = dff.copy() #dff.apply(lambda x: x/x.max(), axis=0).copy()
    normdff['median'] = normdff.median(axis=1)
    vector_dict[i] = list(normdff['median'])



vectordf = pd.DataFrame({'vectorgermline':vector_dict['vectorgermline'],'vectortumor':vector_dict['vectortumor']})
vectordf.to_csv('median_vectors.csv',sep='\t',index=None)
vector_color = {'vectorgermline':'r','vectortumor':'b'}
vector_legend_label = {'vectorgermline':'GERMLINE','vectortumor':'TUMOR'}

for run in list(vectordf.columns):
    vector = list(vectordf[run])#[::20]
    t = np.arange(start=0,stop=len(vector), step = 1)
    plt.plot(t, vector,vector_color[run] , label = vector_legend_label[run], linewidth=1)

plt.legend(loc='upper right',fontsize='medium') #bbox_to_anchor=(1.05, 1),borderaxespad=0.
plt.ylabel('Normalized mtDNA coverage',fontsize=14)
plt.xlabel('mtDNA position',fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.show()

