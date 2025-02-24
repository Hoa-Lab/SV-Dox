from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
from multiprocessing import Pool

#-----------------------------------------------------
fd_out='./out/a01_dd_00_prep'
fd_in='./out/a00_prep_02_qc'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_in).glob('*.h5ad'))

#-----------------------------------------------------
def mainf(fname):
    sample=Path(fname).stem
    ada=sc.read(fname)
    #filter
    ada=ada[ada.obs['qc_cnt']==0, :]
    ada=ada[ada.obs['qc_gene']==0, :]
    ada=ada[ada.obs['qc_mt']==0, :]
    ada.write(f'{fd_out}/{sample}.h5ad') 
    print(ada.shape)
    return

#-----------------------------------------------------
with Pool(10) as p: p.map(mainf, l_fname)
