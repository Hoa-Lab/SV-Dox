from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool

#-----------------------------------------------------
t=0.7
l_sample=['Het0', 'Het1', 'NonF0', 'NonF1', 'Fluc0', 'Fluc1']
d_rename={'raw': 'Raw', 'qc_cnt': 'Filter by Counts', 'qc_gene': 'Filter by Genes', 'qc_mt': 'Filter by Mt%', 'qc_dd': 'Dedoublet'}

fd_out='./out/a02_clean_00_filter'
fd_in='./out/a00_prep_02_qc'
fd_dd='./out/a01_dd_03_r1'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-----------------------------------------------------
def mainf(sample):
    ada=sc.read(f'{fd_in}/{sample}.h5ad')
    ada_dd=sc.read(f'{fd_dd}/{sample}.h5ad')
    #pred
    ada.obs['dd_score']=1
    ada.obs.loc[ada_dd.obs.index, ['dd_score']]=ada_dd.obs['dd_score']
    ada.obs['qc_dd']=(ada.obs['dd_score']>t).astype('int')
    #filter
    l_data=[['raw', ada.shape[0]]]
    for col in ['qc_cnt', 'qc_gene', 'qc_mt', 'qc_dd']:
        ada=ada[ada.obs[col]==0, :].copy()
        l_data.append([col, ada.shape[0]])
    df=pd.DataFrame(l_data, columns=['filter', 'cnt'])
    #save
    ada.write(f'{fd_out}/{sample}.h5ad')
    df.to_csv(f'{fd_out}/{sample}.csv')
    return

def plt_qc(df, f_out, title='Cell Counts After Sequential Filtering', sz=(8,3.5), cmap=['#666666', '#828282', '#999999', '#b0b0b0', '#ff0000']):
    #prep
    if sz is None: sz=(9, 0.1*df.shape[0])
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.barplot(x='sample', y='cnt', data=df, hue='filter', palette=cmap)
    #adjust
    sns.despine()
    ax.set_title(title, x=0.6, fontsize=22, weight='semibold', pad=10)
    plt.xlabel('')
    plt.ylabel('Counts', fontsize=15, labelpad=12, weight='semibold')
    plt.xticks(fontsize=14, weight='semibold')
    plt.legend(loc=(1.01, 0), frameon=False, prop={'weight': 'semibold', 'size':10})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-----------------------------------------------------
#filter
#with Pool(10) as p: p.map(mainf, l_sample)

#plot - prep
l_df=[]
for sample in l_sample:
    df=pd.read_csv(f'{fd_out}/{sample}.csv', index_col=0)
    df['sample']=sample
    l_df.append(df)
df=pd.concat(l_df)
df['filter']=df['filter'].apply(lambda x: d_rename.get(x, x))

#plot
f_out=f'{fd_out}/cnt.png'
plt_qc(df, f_out)
