from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import json

#-----------------------------------------------------
l_update=['Slc4a10+', 'Lrp2+']  #cells need to update 

fd_out='./out/b01_r1_00_prep'
f_ada='./out/b00_r0_01_train/data.h5ad'
f_mark='./out/b00_r0_00_mark/mark.json'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
with open(f_mark, 'r') as f: d_mark=json.load(f)
ada=sc.read(f_ada)

#-----------------------------------------------------
#prep
l_update=l_update+['Unknown']
l_keep=[i for i in ada.obs['predict'].astype('str').unique() if i not in l_update]

#marker blacklist
l_block=[d_mark[i] for i in l_keep]
l_block=sum(l_block, [])
Path(f'{fd_out}/gene_bl.txt').write_text('\n'.join(l_block))

#clean ada
ada=ada[~ada.obs['predict'].isin(l_keep), :].copy()
ada=sc.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
ada.obs=ada.obs.drop(['prelabel', 'predict'], axis=1)

ada.write(f'{fd_out}/data.h5ad')
print(ada.shape)
