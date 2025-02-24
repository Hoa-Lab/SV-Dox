from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from multiprocessing import Pool
import re

#-----------------------------------------------------
l_name=['Het0', 'Het1', 'NonF0', 'NonF1', 'Fluc0', 'Fluc1']

fd_out='./out/a00_prep_00_load'
fd_in='../a00_raw/data'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-----------------------------------------------------
def filter_gene(l_gene, l_ignore=['Ercc-', 'Tg-', 'Tg_', 'n-', '-ps', '-as', '.']):
    #general
    for t in l_ignore: l_gene=[i for i in l_gene if not t in i]  
    l_gene=[i for i in l_gene if not re.match('^Gm[0-9][0-9].*', i)]
    l_gene=[i for i in l_gene if not re.match('^[a-z].*', i)]
    l_gene=[i for i in l_gene if not re.match('^[A-Z][A-Z].*', i)]
    l_gene=[i for i in l_gene if len(i)<15] 
    #other
    l_gene=[i for i in l_gene if not re.match('^Rp[l,s].*', i)]
    l_gene=[i for i in l_gene if not re.match('^Mir[0-9][0-9].*', i)]
    l_gene=[i for i in l_gene if not re.match('^Linc[0-9][0-9].*', i)]
    l_gene=[i for i in l_gene if not (('Rik' in i) & (not i.startswith('Rik')))]
    l_gene=[i for i in l_gene if not 'mt-' in i]
    #sort
    l_gene.sort()
    return l_gene


def load_ada(fd, tag=None, d_extra={}, mt='mt-'):
    #load
    ada=sc.read_10x_mtx(fd, cache=False)
    sc.pp.filter_genes(ada, min_counts=1)
    ada.obs.index=ada.obs.index.map(lambda x: f'{x}_{tag}')
    for i, j in d_extra.items(): ada.obs[i]=j  
    #stats
    ada.var['mt']=ada.var_names.str.startswith(mt)
    sc.pp.calculate_qc_metrics(ada, qc_vars=['mt'], percent_top=None, inplace=True, log1p=False)
    #label gene (0: keep, 1: filter)
    ada.var['filter']=1
    l_gene=filter_gene(ada.var.index.tolist())
    ada.var.loc[ada.var.index.isin(l_gene), ['filter']]=0 
    return ada


def mainf(name):
    fd=f'{fd_in}/SV_{name}'
    d_extra={'sample': name[:-1]}
    ada=load_ada(fd, tag=name, d_extra=d_extra)
    ada.write(f'{fd_out}/{name}.h5ad') 
    print(ada.obs)
    return

#-----------------------------------------------------
with Pool(10) as p: p.map(mainf, l_name)
