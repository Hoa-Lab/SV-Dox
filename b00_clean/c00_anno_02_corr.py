from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns

#-----------------------------------------------------
r0=0.2  #min corr
r1=0.2  #max corr for other cells

fd_out='./out/c00_anno_02_corr'
fd_in='./out/c00_anno_00_train'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
df_pred=pd.read_csv(f'{fd_in}/predict.csv', index_col=0)

#cnt
ada=sc.read(f'{fd_in}/data.h5ad')
ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
df=ada.to_df()
df=df.loc[:, df.sum()>0].copy()

#-----------------------------------------------------
def mainf(cell):
    print(cell)
    #corr
    df_corr=pd.DataFrame(df.corrwith(df_pred[cell], axis=0), columns=['corr'])
    df_corr=df_corr.sort_values('corr', ascending=False)
    df_corr.to_csv(f'{fd_out}/{cell}.csv')
    return

#-----------------------------------------------------
l_cell=[i for i in ada.obs['predict'].unique() if i!='Unknown']

#all genes corr
#for cell in l_cell: mainf(cell)

#-----------------------------------------------------
#cocnat corr
l_df=[]
for cell in l_cell:
    df=pd.read_csv(f'{fd_out}/{cell}.csv', index_col=0) 
    df.columns=[cell]
    l_df.append(df)
    df=pd.concat(l_df, axis=1, join='outer')

#filter by corr
df=df.loc[df.max(axis=1)>r0, :].copy()
df['cell']=df.idxmax(axis=1)

#negative filter
l_df=[]
for cell in df['cell'].unique():
    dfi=df.query('cell==@cell')
    dfi=dfi.sort_values(cell, ascending=False)
    dfi=dfi.drop('cell', axis=1)
    #filter
    df_tmp=dfi.drop(cell, axis=1)
    df_tmp=df_tmp.loc[df_tmp.max(axis=1)<r1, :]
    dfi=dfi.reindex(df_tmp.index)
    #clean
    dfi=dfi.reindex([cell], axis=1)
    dfi.columns=['corr']
    dfi.to_csv(f'{fd_out}/{cell}_filter.csv')
