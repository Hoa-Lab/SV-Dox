from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool
from sklearn.mixture import GaussianMixture

#-----------------------------------------------------
l_sample=['Het0', 'Het1', 'NonF0', 'NonF1', 'Fluc0', 'Fluc1']

fd_out='./out/a02_clean_03_gm'
f_in='./out/a02_clean_01_concat/data.h5ad'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_in)
ada.obs['tmp']=ada.obs.index.str.split('_').str[-1]

#-----------------------------------------------------
def get_gm(df):
    #gm
    l_data=[]
    for i, gene in enumerate(df.columns.tolist()):
        print(f'{i} out of {df.shape[1]}: {gene}')
        ar=df[gene].values.flatten()
        ar=ar[ar!=0]
        #calulate mean
        gm=GaussianMixture(n_components=2, random_state=0).fit(ar.reshape([-1, 1]))
        if not gm.converged_: continue
        #cleam param
        df_param=pd.DataFrame({'mean': gm.means_.flatten(), 'cov': gm.covariances_.flatten(), 'w': gm.weights_.flatten()}) 
        df_param=df_param.sort_values('mean')
        l_data.append([gene]+list(df_param.values.flatten()))
    #make df
    df=pd.DataFrame(l_data, columns=['gene', 'm0', 'cov0', 'w0', 'm1', 'cov1', 'w1'])
    df=df.set_index('gene')
    #clean
    df['ratio_m']=df['m1']/df['m0']
    df=df.sort_values('ratio_m', ascending=False)
    return df

def mainf(sample, min_cell=20, min_cnt=0.5):
    adai=ada[ada.obs['tmp']==sample, :].copy() 
    df=ad.AnnData(adai.raw.X, obs=adai.obs, var=adai.raw.var).to_df()
    #filter genes
    dfi=df.T.copy()
    dfi=dfi.loc[(dfi!=0).sum(axis=1)>min_cell, :].copy()
    dfi=dfi.replace(0, np.NaN)
    dfi=dfi.loc[dfi.mean(axis=1)>min_cnt, :].copy()
    df=df.reindex(dfi.index, axis=1)
    #gm
    df=get_gm(df)
    df.to_csv(f'{fd_out}/{sample}.csv')
    return

def combine_col(d_gm, col, ada=ada, max_na=4): 
    l_df=[]
    for sample, df in d_gm.items():
        dfi=df.reindex([col], axis=1)
        dfi.columns=[sample]
        l_df.append(dfi)
    #concat
    df=pd.concat(l_df, axis=1)
    df=df.loc[df.isna().sum(axis=1)<=max_na, :]
    df[col] =df.mean(axis=1)
    df=df.reindex([col], axis=1)    
    return df

#-----------------------------------------------------
#gm
#with Pool(12) as p: p.map(mainf, l_sample)

#combine - prep
d_gm={}
for sample in l_sample: d_gm[sample]=pd.read_csv(f'{fd_out}/{sample}.csv', index_col=0)

#combine
l_df=[]
for col in ['m0', 'm1', 'cov0', 'cov1', 'w0', 'w1']:
    df=combine_col(d_gm, col)    
    l_df.append(df)
df=pd.concat(l_df, axis=1).dropna()

#clean
df=df.query('m0>0.5')
df['ratio_m']=df['m1']/df['m0']
df=df.sort_values('ratio_m', ascending=False)

df.to_csv(f'{fd_out}/gm.csv')
