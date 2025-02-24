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
n=20

fd_out='./out/d00_plot_00_mark-bar'
fd_in='./out/c00_anno_03_anno/corr'
f_meta='../a00_raw/meta/cell.json'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_in).glob('*.csv'))
with open(f_meta, 'r') as f: d_meta=json.load(f)
d_cmap={d_meta[i]['cell']:d_meta[i]['clr'] for i in list(d_meta)}

#-----------------------------------------------------
def plt_bar(df, f_out, title=None, sz=(4.5, 8), clr='grey', xlim=0.9, adjust=None, xi=0):
    #check
    if df.shape[0]==0: return
    #plot ([left, bottom, width, height])
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.barplot(x='corr', y=df.index, data=df, color=clr, alpha=0.9)
    #adjust
    plt.title(title, fontsize=30, weight='semibold', pad=12, x=0.42)
    plt.xlabel('Marker Score', fontsize=17, weight='semibold')
    plt.yticks(fontsize=17, weight='semibold')
    plt.xlim([0, xlim])
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def mainf(fname):
    name=Path(fname).stem
    cell=name.split('_')[0]
    df=pd.read_csv(fname, index_col=0)
    df=df.iloc[:n, :] 
    #plot
    f_out=f'{fd_out}/{name}.png'
    clr=d_cmap[cell]
    title=cell
    plt_bar(df, f_out, title=title, clr=clr)
    return

#-----------------------------------------------------
with Pool() as p: p.map(mainf, l_fname)
#mainf(l_fname[0])
