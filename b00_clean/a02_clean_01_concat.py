from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy.external as sce

#-----------------------------------------------------
l_sample=['Het0', 'Het1', 'NonF0', 'NonF1', 'Fluc0', 'Fluc1']
c_sample=['#4287f5', '#8ae332', '#f00505']
d_rename={'Het': 'Het NF', 'NonF': 'DE17.5 NF', 'Fluc': 'DE17.5 F'}

fd_out='./out/a02_clean_01_concat'
fd_in='./out/a02_clean_00_filter'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-----------------------------------------------------
def plt_umap(ada, f_out, title=None, cmap='tab20', sz=(2.7, 2.4), s=0.5, obsm='X_umap', col='sample'):
    #prep
    df=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col]=ada.obs[col].astype('str').replace(d_rename)
    df[col]=pd.Categorical(df[col], categories=list(d_rename.values()), ordered=True)    
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0.8)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=14, weight='semibold')
    plt.legend(loc=(0.75, 0.05), frameon=False, prop={'size': 8, 'weight': 'semibold'}, markerscale=1)
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-----------------------------------------------------
'''
#norm
d_ada={}
for sample in l_sample:
    ada=sc.read(f'{fd_in}/{sample}.h5ad')
    #norm
    sc.pp.normalize_total(ada, exclude_highly_expressed=True)
    sc.pp.log1p(ada)
    #filter genes
    ada=ada[:, ada.var['filter']==0].copy()
    d_ada[sample]=ada

#get common genes
l_gene=[set(ada.var.index.tolist()) for ada in d_ada.values()]
l_gene=list(set.intersection(*l_gene))
l_gene.sort()

#concat - prep
l_df, l_obs, l_var=([], [], [])
for sample, ada in d_ada.items():
    ada=ada[:, l_gene].copy()
    #df & obs
    l_df.append(ada.to_df())
    l_obs.append(ada.obs)
    #var
    df_var=ada.var.copy()
    df_var=df_var.reindex(['gene_ids']+['n_counts', 'n_cells_by_counts'], axis=1)
    df_var.columns=[f'{i}_{sample}' for i in df_var.columns]        
    l_var.append(df_var)
    
#concat
df=pd.concat(l_df)
obs=pd.concat(l_obs)
var=pd.concat(l_var, axis=1).reindex(l_gene)
ada=ad.AnnData(df.values, obs=obs, var=var)
ada.obs['sample']=pd.Categorical(ada.obs['sample'], categories=['Het', 'NonF', 'Fluc'], ordered=True)

#HVG & PCA
ada.raw=ada.copy()
sc.experimental.pp.highly_variable_genes(ada, n_top_genes=2000, flavor='pearson_residuals', inplace=True)
ada=ada[:, ada.var.highly_variable]
sc.experimental.pp.normalize_pearson_residuals(ada)
sc.tl.pca(ada, svd_solver='arpack', n_comps=30)

#UMAP
adai=ada.copy()
sc.pp.neighbors(adai)
sc.tl.umap(adai, n_components=2, random_state=42) 
ada.obsm['X_umap']=adai.obsm['X_umap']

#Harmony UMAP
adai=ada.copy()
adai.obs['name']=adai.obs.index.str.split('_').str[-1]
sce.pp.harmony_integrate(adai, 'name')
sc.pp.neighbors(adai, use_rep='X_pca_harmony')
sc.tl.umap(adai, n_components=2, random_state=42)
ada.obsm['X_harmony']=adai.obsm['X_umap']
ada.write(f'{fd_out}/data.h5ad')
'''
#-----------------------------------------------------
ada=sc.read(f'{fd_out}/data.h5ad')
print(ada.shape)

#plot samples
col='sample'
title='Samples'
f_out=f'{fd_out}/sample.png'
plt_umap(ada, f_out, title=title, cmap=c_sample)

f_out=f'{fd_out}/sample_harmony.png'
plt_umap(ada, f_out, title=title, cmap=c_sample, obsm='X_harmony')
