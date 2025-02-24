from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

#-----------------------------------------------------
l_sample=['Het', 'NonF', 'Fluc']
d_sample={'Het': 'Het NF', 'NonF': 'DE17.5 NF', 'Fluc': 'DE17.5 F'}
#l_name=['Het0', 'Het1', 'NonF0', 'NonF1', 'Fluc0', 'Fluc1']

l_param=[
('Genes by Counts', 'n_genes_by_counts', 'Counts'), 
('Total Counts', 'total_counts', 'Counts'), 
('Mt percent', 'pct_counts_mt', 'Percent')
]

fd_out='./out/a00_prep_01_plt-sum'
fd_in='./out/a00_prep_00_load'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-----------------------------------------------------
def plt_sum(df, col, f_out, title=None, sz=(3.5, 4), ylbl=None):
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.boxplot(y=col, x='sample', data=df, color='grey', showfliers=False, linewidth=0.5)
    #adjust
    plt.title(title, fontsize=20, weight='semibold', pad=10)
    plt.xlabel('')
    plt.ylabel(ylbl, fontsize=12, weight='semibold')
    plt.xticks(fontsize=12, weight='semibold', rotation=45)
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-----------------------------------------------------
#prep
l_fname=list(Path(fd_in).glob('*.h5ad'))
l_ada=[sc.read(i) for i in l_fname]
ada=ad.concat(l_ada)

#change sample names
ada.obs['sample']=ada.obs['sample'].astype('str')
ada.obs['sample']=ada.obs['sample'].replace(d_sample)
ada.obs['sample']=pd.Categorical(ada.obs['sample'], categories=list(d_sample.values()), ordered=True)

#plot
for param in l_param: 
    title, col, ylbl=param

    #plot
    f_out=f'{fd_out}/{col}.png'
    plt_sum(ada.obs, col, f_out, title=title, ylbl=ylbl)
