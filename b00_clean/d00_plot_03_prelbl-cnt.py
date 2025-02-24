from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
import matplotlib.pyplot as plt
import seaborn as sns

#---------------------------------------------
fd_out='./out/d00_plot_03_prelbl-cnt'
f_in='./out/c00_anno_03_anno/data.h5ad'
f_meta='../a00_raw/meta/cell.json'

##############################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
with open(f_meta, 'r') as f: d_meta=json.load(f)
d_match={i: d_meta[i]['cell'] for i in d_meta.keys()}

#count cells
ada=sc.read(f_in)
print(ada.shape)  #9268
adai=ada[ada.obs['prelabel']!='TBD', :]
print(adai.shape) #715

#-----------------------------------------------------
def plt_bar(df, f_out, title=None, cmap=None):
    #plot
    fig, ax=plt.subplots(figsize=(2.5, 3.2))
    ax=sns.barplot(data=df, x='cnt', y=df.index, palette=cmap)
    #adjust
    plt.title(title, fontsize=14, weight='semibold', pad=8, x=0.32)
    plt.xlabel('Counts', fontsize=8, weight='semibold', labelpad=4)
    plt.ylabel('')
    plt.yticks(fontsize=8, weight='semibold')
    plt.xticks(fontsize=5)
    plt.xlim([0, 400])
    #text
    for i in range(df.shape[0]):
        c=df['cnt'].tolist()[i]
        txt=str(c)
        x=c+20
        y=i+0.1
        plt.text(x, y, txt, fontsize=5, weight='semibold', ma='center')
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def plt_umap(ada, f_out, title=None, cmap='tab20', sz=(8,6), s=3, obsm='X_harmony', col='prelabel'):
    #prep
    df=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col]=ada.obs[col]
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0.9)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=27, weight='semibold')
    plt.legend(loc=(1, 0.05), ncol=1, frameon=False, prop={'size': 12, 'weight': 'semibold'})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#---------------------------------------------
#prep
df=ada.obs.reindex(['prelabel'], axis=1)
df['cnt']=1
df=df.groupby('prelabel').count()

#sort
l_cell=list(d_meta.keys())[:-1]
l_cat=[d_meta[i]['cell'] for i in l_cell]
df=df.reindex(l_cell)
df.index=l_cat

#plot bar
f_out=f'{fd_out}/prelabel.png'
title='Prelabeled Cells'
cmap=[d_meta[i]['clr'] for i in l_cell]
plt_bar(df, f_out, cmap=cmap, title=title)

#---------------------------------------------
#prep umap
#ada.obs['tmp']=ada.obs['prelabel'].astype('str').apply(lambda x: d_match.get(x, 'TBD'))
#ada.obs['tmp']=pd.Categorical(ada.obs['tmp'], categories=l_cat+['TBD'], ordered=True)

#plot umap
#f_out=f'{fd_out}/prelabel_umap.png'
#cmap=[d_meta[i]['clr'] for i in l_cell]+['#d4d6d9']
#plt_umap(ada, f_out, title='Prelabeled Cells', cmap=cmap, col='tmp')
