from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import json
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns

#-----------------------------------------------------
l_gene=['Esrrb', 'Tyr', 'Cldn11', 'Anxa1', 'Slc26a4', 'Epyc', 'Lgr5', 'Coch', 'Otos', 'Runx2', 'Meg3', 'Slc26a7', 'Stac', 'Rgs5', 'Hba-a1', 'Hba-a2', 'Pdzrn4', 'Nrxn3', 'Ptprz1', 'Fgf12', 'Nkain2', 'Pgm5', 'Kcnq1ot1']
l_gene=['Pid1', 'Pde3a']

fd_out='./out/a02_clean_02_exp'
f_in='./out/a02_clean_01_concat/data.h5ad'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_in)

#-----------------------------------------------------
def plt_exp(gene, ada, f_out, sz=(7,6), s=3.5, cmap='BuPu', vmin=-0.6, vmax=None, obsm='X_harmony', raw=True):
    #prep 
    df_xy=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var).to_df()
    df=df.merge(df_xy, left_index=True, right_index=True)
    #norm cmap
    if vmax is None: vmax=df[gene].max()
    norm=plt.Normalize(vmin, vmax)
    l_tick=np.arange(-0.5, vmax, 0.5).tolist()
    #plot  (left, bottom, width, height)
    fig=plt.figure(figsize=sz)
    ax=fig.add_axes([0, 0, 0.85, 0.9])
    ax=sns.scatterplot(x='x', y='y', data=df, hue=gene, palette=cmap, s=s, linewidth=0, alpha=0.9, hue_norm=norm)
    #adjust
    plt.axis('off')
    ax.set_title(gene, fontsize=34, weight='semibold')
    ax.get_legend().remove()
    #color bar (https://stackoverflow.com/questions/62884183/trying-to-add-a-colorbar-to-a-seaborn-scatterplot)
    ax_cbar=fig.add_axes([0.85, 0.3, 0.02, 0.4])
    sm=plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb=fig.colorbar(sm, cax=ax_cbar, ticks=l_tick)
    cb.ax.tick_params(labelsize=14)
    #save
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

def mainf(gene):
    f_out=f'{fd_out}/{gene}.png'
    try: plt_exp(gene, ada, f_out)
    except Exception: print(f'{gene} error...')
    return

#-----------------------------------------------------
with Pool(10) as p: p.map(mainf, l_gene)

#---------------------------------------------
gene='Anxa1'
f_out=f'{fd_out}/{gene}.png'
#plt_exp(gene, ada, f_out, obsm='X_harmony', vmin=-0.5, vmax=2.4)
