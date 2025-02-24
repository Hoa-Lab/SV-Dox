from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from matplotlib.patches import Patch
import json

#-----------------------------------------------------
l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle']
c_cell=['#F8766D', '#00BA38', '#619CFF', '#CD950C']
l_sample=['Het', 'NonF', 'Fluc']
c_sample=['#4287f5', '#8ae332', '#f00505']
n=10

fd_out=f'./out/d00_plot_01_mark-hm'
fd_in='./out/c00_anno_03_anno/corr'
f_ada='./out/c00_anno_03_anno/data.h5ad'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#count
ada=sc.read(f_ada)
ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
df_cnt=ada.to_df()
df_cnt['cell']=ada.obs['cell'].astype('str')
df_cnt['sample']=ada.obs['sample']

#-----------------------------------------------------
def plt_hm(df, df_sample, f_out, title=None, sz=(10, 9)):
    fig=plt.figure(figsize=sz) 
    #plot top cell bar [left, bottom, width, height]
    ax_top=fig.add_axes([0.1, 0.87, 0.66, 0.05])
    ax_top=sns.heatmap(df_sample, cmap=c_cell+c_sample, cbar=False)
    ax_top.axes.get_xaxis().set_visible(False)
    ax_top.set_title(title, fontsize=27, pad=10, weight='semibold')
    plt.yticks(fontsize=12, rotation=0, weight='semibold', ticks=[0.5, 1.5], labels=['Cell', 'Sample']) 
    #plot data
    ax=fig.add_axes([0.1, 0.02, 0.66, 0.835])
    ax=sns.heatmap(df, cmap='coolwarm', cbar=False, vmax=2, vmin=-2, yticklabels=True)
    ax.axes.get_xaxis().set_visible(False)
    plt.yticks(fontsize=12.5, rotation=0, weight='semibold')
    #legend cell
    l_legend=[]
    for i in range(len(l_cell)):
        clr, lbl=c_cell[i], l_cell[i]
        p=Patch(facecolor=clr, edgecolor=clr, label=lbl, lw=3)
        l_legend.append(p)
    ax_legend=fig.add_axes([0.935, 0.625, 0.04, 0.45])
    ax_legend.set_axis_off()
    ax_legend.legend(handles=l_legend, loc='right', title='Cell', prop={'size': 12, 'weight': 'semibold'}, title_fontproperties={'size':12, 'weight': 'semibold'})
    #legend sample
    l_legend=[]
    for i in range(len(l_sample)):
        clr, lbl=c_sample[i], l_sample[i]
        p=Patch(facecolor=clr, edgecolor=clr, label=lbl, lw=3)
        l_legend.append(p)
    ax_legend=fig.add_axes([0.861, 0.46, 0.04, 0.45])
    ax_legend.set_axis_off()
    ax_legend.legend(handles=l_legend, loc='right', title='Sample', prop={'size': 12, 'weight': 'semibold'}, title_fontproperties={'size':12, 'weight': 'semibold'})
    #cbar
    ax_cbar=fig.add_axes([0.81, 0.04, 0.02, 0.45])
    norm=plt.Normalize(-2, 2)
    sm=plt.cm.ScalarMappable(cmap='coolwarm', norm=norm)
    cb=fig.colorbar(sm, cax=ax_cbar, ticks=[-2, -1, 0, 1, 2])
    cb.set_label('Normalized Counts (z-score)', rotation=270, weight='semibold', fontsize=12.5, labelpad=15)
    cb.ax.yaxis.tick_left()
    cb.ax.tick_params(labelsize=10) 


    #save
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-----------------------------------------------------
#gene
l_gene=[]
for cell in l_cell:
    df=pd.read_csv(f'{fd_in}/{cell}.csv', index_col=0)
    l_gene.extend(df.index.tolist()[:n]) 

#cnt
dfc=df_cnt.reindex(l_gene+['cell', 'sample'], axis=1)
dfc=dfc.loc[dfc['cell'].isin(l_cell), :].copy()
dfc['cell']=pd.Categorical(dfc['cell'], categories=l_cell, ordered=True)
dfc=dfc.sort_values(['cell', 'sample'])

#prep
df_sample=dfc.reindex(['cell', 'sample'], axis=1).T
df_sample=df_sample.replace(l_cell+l_sample, range(len(l_cell)+len(l_sample)))
dfc=dfc.drop(['cell', 'sample'], axis=1).T
dfc=dfc.apply(zscore, axis=1)

#plot
f_out=f'{fd_out}/mark.png'
title='Top Markers in Major SV Cells'
plt_hm(dfc, df_sample, f_out, title=title)

