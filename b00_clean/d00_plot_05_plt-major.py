from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import json

#---------------------------------------------
l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle']

fd_out='./out/d00_plot_05_plt-major'
f_in='./out/c00_anno_03_anno/data.h5ad'
f_meta='../a00_raw/meta/cell.json'

##############################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
with open(f_meta, 'r') as f: d_meta=json.load(f)
ada=sc.read(f_in)

#---------------------------------------------
def plt_umap(ada, f_out, title=None, sz=(7.5,6), s=3.5, obsm='X_harmony', col='tmp'):
    #prep
    df=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col]=ada.obs[col]
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0.8)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=27, weight='semibold')
    plt.legend(loc=(0.95, 0.3), frameon=False, prop={'size': 14, 'weight': 'semibold'})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#---------------------------------------------
#prep
ada.obs['tmp']=ada.obs['cell'].astype('str')
ada.obs['tmp']=ada.obs['tmp'].apply(lambda x: x if (x in l_cell) else 'Others')
ada.obs['tmp']=pd.Categorical(ada.obs['tmp'], categories=l_cell+['Others'], ordered=True)

#plot
d_cmap={d_meta[i]['cell']:d_meta[i]['clr'] for i in d_meta.keys()}
cmap=[d_cmap[i] for i in l_cell]+['#d4d6d9']

f_out=f'{fd_out}/cell.png'
title='Major SV Cells'
plt_umap(ada, f_out, title=title)
