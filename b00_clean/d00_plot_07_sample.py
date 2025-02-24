from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import json
from multiprocessing import Pool

#-----------------------------------------------------
l_sample=['Het', 'NonF', 'Fluc']
d_sample={'Het': 'Het NF', 'NonF': 'DE17.5 NF', 'Fluc': 'DE17.5 F'}

fd_out='./out/d00_plot_07_sample'
f_in='./out/a02_clean_01_concat/data.h5ad'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_in)

#-----------------------------------------------------
def plt_umap(ada, f_out, title=None, cmap='tab20', sz=(7,6), s=3.5, obsm='X_harmony', col='sample'):
    #prep
    df=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col]=ada.obs[col]
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0.5)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=27, weight='semibold')
    plt.legend(loc=(0.9, 0.05), frameon=False, prop={'size': 16, 'weight': 'semibold'})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-----------------------------------------------------
for sample in l_sample:
    ada.obs['tmp']=ada.obs['sample'].replace([i for i in l_sample if i!=sample], 'Others')
    ada.obs['tmp']=ada.obs['tmp'].replace(sample, d_sample[sample])
    ada.obs['tmp']=pd.Categorical(ada.obs['tmp'], categories=[d_sample[sample], 'Others'], ordered=True)
    #plot
    f_out=f'{fd_out}/{sample}.png'
    plt_umap(ada, f_out, col='tmp', cmap=['#ff0000', '#999999'], title=d_sample[sample])
    

