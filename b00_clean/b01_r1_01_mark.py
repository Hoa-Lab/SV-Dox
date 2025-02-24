from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from scipy.stats import zscore

#-----------------------------------------------------
max_gene=400
min_mark=3
n_prelbl=2  #gene number for prelabel
z_mark=2
z_other=2
min_cell=5

min_corr=0.3
min_ratio_m=2
min_w=0
l_extra=[]

fd_out='./out/b01_r1_01_mark'
fd_in='./out/b01_r1_00_prep'
f_gm='./out/a02_clean_03_gm/gm.csv'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#gm
l_block=Path(f'{fd_in}/gene_bl.txt').read_text().split('\n')
df=pd.read_csv(f_gm, index_col=0)
df=df.loc[~df.index.isin(l_block), :].copy()

#ada
ada=sc.read(f'{fd_in}/data.h5ad')
ada_save=ada.copy()
df_cnt=ada.to_df()

#-----------------------------------------------------
def gm_group(l_gene, df, min_corr=min_corr):
    #prep
    l_save=l_gene.copy() 
    l_grp=[]
    #get genes
    while len(l_gene)>0:
        gene=l_gene[0]  
        l_grp.append(gene)
        #remove low corr genes (at least 1 gene remains, since corr with self is 1)
        l_tmp=l_gene.copy()
        for g in l_gene: 
            if (df[gene].corr(df[g]))<min_corr: l_tmp=[i for i in l_tmp if i!=g] 
        #update
        l_gene=l_tmp[1:].copy()
    #update remaining gene list
    l_gene=[i for i in l_save if i not in l_grp]
    return l_grp, l_gene

#-----------------------------------------------------
#filter gm
df=df.query('ratio_m>@min_ratio_m')
df=df.query('w0>@min_w')
df=df.query('w1>@min_w')
df=df.sort_values('ratio_m', ascending=False)
df.to_csv(f'{fd_out}/gm_filtered.csv')
print(f'filtered gm: {df.shape}')

#group markers
d_mark={}
l_gene=l_extra+df.index.tolist()[:max_gene]
l_gene=list(dict.fromkeys(l_gene))   #in case duplicated genes

while len(l_gene)>=min_mark:
    l_grp, l_gene=gm_group(l_gene, df_cnt)
    if len(l_grp)>=min_mark: d_mark[f'{l_grp[0]}+']=l_grp
with open(f'{fd_out}/mark.json', 'w') as f: json.dump(d_mark, f)
print(list(d_mark.keys()))

#prelabel-prep
l_sample=ada.obs['sample'].unique().tolist()
l_all=[d_mark[i][:n_prelbl] for i in list(d_mark)]
l_all=sum(l_all, [])

#prelabel
l_df=[]
for sample in l_sample: #prelabel is applied on each sample
    #prep
    df=ada[ada.obs['sample']==sample, :].to_df()
    df=df.reindex(l_all, axis=1)
    df=df.apply(zscore)
    df['prelabel']='TBD'

    #label
    for cell in d_mark.keys():
        l_mark=d_mark[cell][:n_prelbl]
        l_other=[i for i in l_all if i not in l_mark]
        dfi=df.copy()
        #filter
        for g in l_mark: dfi=dfi.loc[dfi[g]>z_mark, :]
        for g in l_other: dfi=dfi.loc[dfi[g]<z_other, :]
        df.loc[dfi.index, ['prelabel']]=cell
    df=df.reindex(['prelabel'], axis=1)
    l_df.append(df)

#prelabel-clean (find low cell number clus)
df=pd.concat(l_df)
df_save=df.copy()

df['cnt']=1
df=df.groupby('prelabel').count()
df=df.query('cnt>@min_cell')
df=df.sort_values('cnt', ascending=False)
print(df)
print(df.shape)

#prelabel-update ada
l_cell=df.index.tolist()
df_save['prelabel']=df_save['prelabel'].apply(lambda x: x if x in l_cell else 'TBD')
ada_save.obs=ada_save.obs.merge(df_save, left_index=True, right_index=True)
ada_save.write(f'{fd_out}/data.h5ad')
