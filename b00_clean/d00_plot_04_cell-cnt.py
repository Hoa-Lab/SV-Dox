from pathlib import Path
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import json
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

#-----------------------------------------------------
l_sample=['Het', 'NonF', 'Fluc']
c_sample=['#4287f5', '#8ae332', '#f00505']
d_rename={'Het': 'Het NF', 'NonF': 'DE17.5 NF', 'Fluc': 'DE17.5 F'}

fd_out='./out/d00_plot_04_cell-cnt'
f_in='./out/c00_anno_03_anno/data.h5ad'
f_cell='../a00_raw/meta/cell.json'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
with open(f_cell, 'r') as f: d_cell=json.load(f)
ada=sc.read(f_in)

#-----------------------------------------------------
def plt_bar(df, f_out, title=None, cmap=None):
    #plot
    fig, ax=plt.subplots(figsize=(2.5, 3.4))
    ax=sns.barplot(data=df, x='cnt', y=df.index, palette=cmap)
    #adjust
    plt.title(title, fontsize=14, weight='semibold', pad=8)
    plt.xlabel('Counts', fontsize=8, weight='semibold', labelpad=4)
    plt.ylabel('')
    plt.xticks([0, 1000, 2000], fontsize=6)
    plt.yticks(fontsize=8, weight='semibold')
    plt.xlim([0, 2200])
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

def plt_cumsum(df, f_out, sz=(2.8, 3.2), l_sample=None, cmap=None):
    df=df.T
    df=df*100
    #plot  [left, bottom, width, height]
    fig=plt.figure(figsize=sz)
    ax=fig.add_axes([0.34, 0.14, 0.43, 0.76])
    for sample, clr in zip(l_sample[::-1], cmap[::-1]):
        ax=sns.barplot(data=df, y=df.index, x=sample, color=clr)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    #adjust
    plt.title('Sample Distributions', fontsize=14, pad=8, weight='semibold', x=0.4)
    plt.ylabel('')
    plt.xlabel('Percent', fontsize=8, weight='semibold', labelpad=4)
    plt.xticks([0, 50, 100], fontsize=6)
    plt.yticks(fontsize=8, weight='semibold')
    #cell legend (https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html)
    l_legend=[]
    for sample, clr in zip(l_sample, cmap):
        sample=sample if (sample!='Ctrl') else 'Het'
        sample=d_rename[sample]
        p=Patch(facecolor=clr, edgecolor=clr, label=sample, lw=3)
        l_legend.append(p)
    ax_l=fig.add_axes([0.965, 0.62, 0.04, 0.45])
    ax_l.set_axis_off()
    ax_l.legend(handles=l_legend, loc='right', prop={'size': 5, 'weight': 'semibold'})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return


#-----------------------------------------------------
#prep
df=ada.obs.reindex(['cell'], axis=1)
df['cnt']=1
df=df.groupby('cell').count()
df=df.iloc[:-1, :]  #not show Unknown
df.index=df.index.astype('str')

#plot
f_out=f'{fd_out}/cell.png'
title='Cell Counts'
d_cmap={d_cell[i]['cell']: d_cell[i]['clr'] for i in d_cell.keys()}
cmap=[d_cmap[i] for i in df.index]
plt_bar(df, f_out, title=title, cmap=cmap)

#-----------------------------------------------------
#count percentage

#prep
df=ada.obs.reindex(['sample', 'cell'], axis=1)
l_cell=df['cell'].cat.categories.tolist()[:-1]
l_df=[]
for cell in l_cell:
    dfi=df.query('cell==@cell').copy()
    dfi['sample']=pd.Categorical(dfi['sample'], categories=l_sample, ordered=True)
    dfi=dfi.groupby('sample').count()
    dfi['cell']=dfi['cell']/dfi['cell'].sum()
    dfi.columns=[cell]
    l_df.append(dfi)
df=pd.concat(l_df, axis=1)
df.to_csv(f'{fd_out}/pct.csv')

#plot
df=df.cumsum()
f_out=f'{fd_out}/pct.png'
plt_cumsum(df, f_out, l_sample=l_sample, cmap=c_sample)
