from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns

#-----------------------------------------------------
fd_out='./out/c00_anno_01_clus'
f_in='./out/c00_anno_00_train/data.h5ad'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_in)

#-----------------------------------------------------
def plt_clus(cell, ada, f_out, title=None, cmap=['#ff0000', '#888888'], obsm='X_harmony', sz=(7,6), s=3.5):
    #prep
    adai=ada.copy()
    adai.obs['clus']=adai.obs['predict'].apply(lambda x: cell if x==cell else 'Others')
    #prep
    df=pd.DataFrame(adai.obsm[obsm], columns=['x', 'y'], index=adai.obs.index)
    df['clus']=pd.Categorical(adai.obs['clus'], categories=[cell, 'Others'], ordered=True)
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue='clus', palette=cmap, s=s, linewidth=0, alpha=0.8)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=24, weight='semibold')
    plt.legend(loc=(0.9, 0.05), frameon=False, prop={'size': 16, 'weight': 'semibold'})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def mainf(cell):
    f_out=f'{fd_out}/{cell}.png' 
    plt_clus(cell, ada, f_out, title=cell)
    return

#-----------------------------------------------------
l_cell=ada.obs['predict'].unique().tolist()
#with Pool(4) as p: p.map(mainf, l_cell)
for cell in l_cell: mainf(cell)
