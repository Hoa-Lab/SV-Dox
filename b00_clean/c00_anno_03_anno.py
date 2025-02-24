from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import json

#-----------------------------------------------------
fd_out='./out/c00_anno_03_anno'
f_ada='./out/c00_anno_00_train/data.h5ad'
fd_corr='./out/c00_anno_02_corr'
f_cell='../a00_raw/meta/cell.json'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
with open(f_cell, 'r') as f: d_cell=json.load(f)
ada=sc.read(f_ada)

#-----------------------------------------------------
def plt_umap(ada, f_out, title='Cell Annotation', cmap='tab20', sz=(2.8, 2.6), s=0.8, obsm='X_harmony', col='cell'):
    #prep
    df=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col]=ada.obs[col]
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0.8)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=14, weight='semibold')
    ax.get_legend().remove()
    #plt.legend(loc=(1.05, 0.05), ncol=1, markerscale=0.7, frameon=False, prop={'size': 3, 'weight': 'semibold'})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def plt_legend(ada, f_out, title='Cell Annotation', cmap='tab20', sz=(3, 2.6), s=0.5, obsm='X_harmony', col='cell'):
    #prep
    df=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col]=ada.obs[col]
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=14, weight='semibold')
    plt.legend(loc=(-0.1, 0), ncol=3, markerscale=0.7, frameon=False, prop={'size': 6, 'weight': 'semibold'})
    #save
    #plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def plt_clus(cell, ada, f_out, title=None, cmap=['#ff0000', '#888888'], obsm='X_harmony', sz=(7,6), s=3.5):
    #prep
    adai=ada.copy()
    adai.obs['clus']=adai.obs['cell'].apply(lambda x: cell if x==cell else 'Others')
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

#-----------------------------------------------------
'''
#rename cell
l_cell=[d_cell[i]['cell'] for i in d_cell.keys()]
ada.obs['cell']=ada.obs['predict'].apply(lambda x: d_cell[x]['cell'])
ada.obs['cell']=pd.Categorical(ada.obs['cell'], categories=l_cell, ordered=True)
ada.write(f'{fd_out}/data.h5ad')
'''

#plot cells
ada=sc.read(f'{fd_out}/data.h5ad')
f_out=f'{fd_out}/cell_harmony.png'
cmap=[d_cell[i]['clr'] for i in d_cell.keys()]
plt_umap(ada, f_out, cmap=cmap)

f_out=f'{fd_out}/cell.png'
#plt_umap(ada, f_out, cmap=cmap, obsm='X_umap')

f_out=f'{fd_out}/legend.png'
#plt_legend(ada, f_out, cmap=cmap, obsm='X_umap')
'''
#-----------------------------------------------------
#rename corr
Path(f'{fd_out}/corr').mkdir(exist_ok=True, parents=True)
for clus in d_cell.keys():
    if clus=='Unknown': continue
    cell=d_cell[clus]['cell']
    #all    
    df=pd.read_csv(f'{fd_corr}/{clus}.csv', index_col=0)  
    df.to_csv(f'{fd_out}/corr/{cell}.csv')
    #filtered
    df=pd.read_csv(f'{fd_corr}/{clus}_filter.csv', index_col=0)  
    df.to_csv(f'{fd_out}/corr/{cell}_filter.csv')

#---------------------------------------------
#umap
ada=sc.read(f'{fd_out}/data.h5ad')
Path(f'{fd_out}/clus').mkdir(exist_ok=True, parents=True)
l_cell=ada.obs['cell'].unique().tolist()
for cell in l_cell:
    f_out=f'{fd_out}/clus/{cell}.png' 
    plt_clus(cell, ada, f_out, title=cell)
'''
