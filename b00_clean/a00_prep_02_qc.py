from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool

#-----------------------------------------------------
fd_out='./out/a00_prep_02_qc'
fd_in='./out/a00_prep_00_load'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_in).glob('*.h5ad'))

#-----------------------------------------------------
def qc_ada(ada, q0_cnt=0.01, q1_cnt=0.99, q0_gene=0.01, q1_gene=0.99, min_mt=10, max_cnt=10000, max_gene=5000):
    #prep
    min_cnt=int(ada.obs['total_counts'].quantile(q0_cnt))
    max_cnt=min(int(ada.obs['total_counts'].quantile(q1_cnt)), max_cnt) 
    min_gene=int(ada.obs['n_genes_by_counts'].quantile(q0_gene))
    max_gene=min(int(ada.obs['n_genes_by_counts'].quantile(q1_gene)), max_gene)
    #QC (0: keep, 1: filter)
    ada.obs['qc_cnt']=0
    ada.obs['qc_gene']=0
    ada.obs['qc_mt']=0
    ada.obs.loc[(ada.obs['total_counts']<min_cnt) | (ada.obs['total_counts']>max_cnt), ['qc_cnt']]=1
    ada.obs.loc[(ada.obs['n_genes_by_counts']<min_gene) | (ada.obs['n_genes_by_counts']>max_gene), ['qc_gene']]=1
    ada.obs.loc[ada.obs['pct_counts_mt']>min_mt, ['qc_mt']]=1
    return ada

def mainf(fname):
    sample=Path(fname).stem
    ada=sc.read(fname)
    ada=qc_ada(ada)
    ada.write(f'{fd_out}/{sample}.h5ad')
    return

#-----------------------------------------------------
with Pool(10) as p: p.map(mainf, l_fname)

'''
#check total cell number before dd
l_fname=list(Path(fd_out).glob('*.h5ad'))
l=[]
for fname in l_fname:
    ada=sc.read(fname)
    l.append(ada.shape[0])
print(sum(l)) #12400
'''
