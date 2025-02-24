from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool

#-----------------------------------------------------
fd_out='./out/a01_dd_04_r1-plt'
fd_in='./out/a01_dd_03_r1'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_in).glob('*.h5ad'))

#-----------------------------------------------------
def plt_pred(ada, f_out, title=None, sz=(5,3), n_bin=100):
    #plot
    fig, ax=plt.subplots(figsize=sz)
    sns.despine()
    ax=sns.histplot(ada.obs, x='dd_score', bins=n_bin, color='grey')
    #adjust
    plt.title(title, fontsize=18, weight='semibold', pad=10)
    plt.xlabel('')
    plt.ylabel('Cell Counts', fontsize=13, weight='semibold')
    plt.xlim([0,1])
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def mainf(fname):
    sample=Path(fname).stem
    ada=sc.read(fname)
    #plot
    f_out=f'{fd_out}/{sample}.png'
    plt_pred(ada, f_out, title=sample)
    return


#-----------------------------------------------------
with Pool(10) as p: p.map(mainf, l_fname)
