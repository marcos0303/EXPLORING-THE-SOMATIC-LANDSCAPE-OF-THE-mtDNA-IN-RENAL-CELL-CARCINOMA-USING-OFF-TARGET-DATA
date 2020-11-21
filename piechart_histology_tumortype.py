import matplotlib.pyplot as plt
import matplotlib as mt
df = pd.read_csv('RCC_samples_mtDNA_db_PAIRED_available_20201019.csv',sep='\t')


group_names = list(df.groupby('histology_simplified').count()['sample'].index)
group_size=list(df.groupby('histology_simplified').count()['sample'].values)
histologylabels = [j +  ' \n ' + k + '%' for j,k in zip(group_names,[str(i) for i in list(np.round(np.array(group_size) / np.sum(group_size) *100, decimals = 1))])]
subgroup_names=['.'.join(i) for i in list(df.groupby(['histology_simplified','type_simplified']).count()['sample'].index)]
subgroup_names=[i[1] for i in list(df.groupby(['histology_simplified','type_simplified']).count()['sample'].index)]
subgroup_size=list(df.groupby(['histology_simplified','type_simplified']).count()['sample'].values)

# Create colors
a, b, c, d=[plt.cm.Reds, plt.cm.Blues, plt.cm.Greens, plt.cm.Purples]
# font
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 12}

mt.rc('font', **font)
titlefontdict={'fontsize': 16,
 'fontweight' : 'bold',
 'verticalalignment': 'baseline',
 'horizontalalignment': 'center'}


# First Ring (outside)
fig, ax = plt.subplots()
ax.axis('equal')
mypie, _ = ax.pie(group_size, radius=1.3, labels=histologylabels, colors=[a(0.8), b(0.8), c(0.8), d(0.8)] )
plt.setp( mypie, width=0.3, edgecolor='white')

# Second Ring (Inside)
mypie2, _ = ax.pie(subgroup_size, radius=1.3-0.3, labeldistance=0.7, colors=[a(0.8), b(0.2), b(0.5), b(0.1), c(0.2), c(0.5), d(0.5), d(0.2)]) #, labels=subgroup_names
plt.setp( mypie2, width=0.4, edgecolor='white')
plt.margins(0,0)
plt.title('RCC histology subtypes (n=142)', fontdict=titlefontdict, loc='center', pad=25)
# show it
plt.show()
