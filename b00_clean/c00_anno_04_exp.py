from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool

#-----------------------------------------------------
n=5

fd_out='./out/c00_anno_04_exp'
fd_in='./out/c00_anno_03_anno/corr'
f_ada='./out/c00_anno_00_train/data.h5ad'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_in).glob('*_filter.csv'))
ada=sc.read(f_ada)

#-----------------------------------------------------
def plt_exp(gene, ada, f_out, sz=(7,6), s=3.5, cmap='BuPu', vmin=-0.6, vmax=None, obsm='X_harmony', raw=True):
    #prep 
    df_xy=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    if raw: df=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var).to_df()
    else: df=ada.to_df()
    df=df.merge(df_xy, left_index=True, right_index=True)
    #norm cmap
    if vmax is None: vmax=df[gene].quantile(1)
    norm=plt.Normalize(vmin, vmax)
    l_tick=np.arange(-0.5, vmax, 0.5).tolist()
    #plot  (left, bottom, width, height)
    fig=plt.figure(figsize=sz)
    ax=fig.add_axes([0, 0, 0.85, 0.9])
    ax=sns.scatterplot(x='x', y='y', data=df, hue=gene, palette=cmap, s=s, linewidth=0, alpha=0.9, hue_norm=norm)
    #adjust
    plt.axis('off')
    ax.set_title(gene, fontsize=30, weight='semibold')
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

def mainf(fname):
    cell=Path(fname).stem.split('_')[0]
    Path(f'{fd_out}/{cell}').mkdir(exist_ok=True, parents=True)
    #get genes
    df=pd.read_csv(fname, index_col=0)
    for gene in df.index.tolist()[:n]:
        f_out=f'{fd_out}/{cell}/{gene}.png'
        plt_exp(gene, ada, f_out)
    return

#-----------------------------------------------------
#for fname in l_fname: mainf(fname)
with Pool(20) as p: p.map(mainf, l_fname)
