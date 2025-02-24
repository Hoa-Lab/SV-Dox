from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns

#---------------------------------------------
l_gene=['CACNA1D', 'CACNA2D1', 'CACNB2', 'CACNG4', 'CACNA1C',
'CACNA1S',
'CACNA1A',
'CACNB1', 
'CACNB3', 
'CACNB4',
'CACNA1G',
'CACNA1H',
'CACNA1I',
'CACNG1-8',
'CACNG2', 'CACNG5',
'CACNG7'
]

#l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root-1', 'Root-2', 'Fibrocyte-1', 'Fibrocyte-2', 'Reissner', 'Pericyte', 'Schwann', 'Endothelial', 'Neuron', 'BMC', 'MMC', 'SC', 'HC', 'RBC']

fd_out='./out/d00_plot_02_mark-violin'
f_in='./out/c00_anno_03_anno/data.h5ad'
f_meta='../a00_raw/meta/cell.json'

##############################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_in)
l_gene=[i.capitalize() for i in l_gene]

#cmap
with open(f_meta, 'r') as f: d_meta=json.load(f)
l_lbl=[i for i in d_meta.keys() if i!='Unknown']
l_cell=[d_meta[i]['cell'] for i in l_lbl]
cmap=[d_meta[i]['clr'] for i in l_lbl]

#-----------------------------------------------------
def plt_violin(ada=None, x=None, y=None, f_out=None, title=None, key='cell', sz=(8,2), raw=True, ignore_zero=False, cmap='tab20', ylabel='Normalized Counts', r=45, hide_x=True):
    #prep
    if raw: ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
    df=ada.to_df().reindex([y], axis=1)
    df[key]=ada.obs[key]
    df=df.loc[df[key].isin(x), :].copy()
    df[key]=pd.Categorical(df[key], categories=x, ordered=True)
    #plot
    fig, ax=plt.subplots(figsize=sz)
    sns.despine()
    ax=sns.violinplot(data=df, x=key, y=y, linewidth=0.5, scale='width', palette=cmap)
    #adjust
    #plt.title(title, fontsize=25, weight='semibold')
    plt.xlabel('')
    plt.ylabel(ylabel, fontsize=17, weight='semibold')
    plt.xticks(fontsize=14, weight='semibold', rotation=r, ha='right', rotation_mode='anchor')
    if hide_x: ax.set_xticks([])
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#---------------------------------------------
for gene in l_gene[:-1]:
    f_out=f'{fd_out}/{gene}.png'
    plt_violin(ada=ada, x=l_cell, y=gene, f_out=f_out, title=gene, key='cell', cmap=cmap, ylabel=gene, hide_x=False, sz=(8, 5))

'''
#last
gene=l_gene[-1]
f_out=f'{fd_out}/{gene}.png'
plt_violin(ada=ada, x=l_cell, y=gene, f_out=f_out, title=gene, key='cell', cmap=cmap, hide_x=False, sz=(8,3), ylabel=gene)
'''
