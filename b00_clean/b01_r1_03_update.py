from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import json

#-----------------------------------------------------
l_update=['Slc4a10+', 'Lrp2+'] 

fd_out='./out/b01_r1_03_update'
f_ada='./out/b01_r1_02_train/data.h5ad'
f_mark='./out/b01_r1_01_mark/mark.json'

f_ada0='./out/b00_r0_01_train/data.h5ad'
f_mark0='./out/b00_r0_00_mark/mark.json'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#mark
with open(f_mark, 'r') as f: d_mark=json.load(f)
with open(f_mark0, 'r') as f: d_mark0=json.load(f)

#ada
ada=sc.read(f_ada)
ada0=sc.read(f_ada0)
ada0.obs['predict']=ada0.obs['predict'].astype('str')
ada0.obs['prelabel']=ada0.obs['prelabel'].astype('str')

#-----------------------------------------------------
def plt_dr(ada, f_out, title=None, cmap='tab20', sz=(8,6), s=4, obsm='X_harmony', col='predict'):
    #prep
    df=pd.DataFrame(ada.obsm[obsm], columns=['x', 'y'], index=ada.obs.index)
    df[col]=ada.obs[col]
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.scatterplot(x='x', y='y', data=df, hue=col, palette=cmap, s=s, linewidth=0, alpha=0.8)
    #adjust
    plt.axis('off')
    ax.set_title(title, fontsize=24, weight='semibold')
    plt.legend(loc=(1.0, 0.05), frameon=False, markerscale=1, prop={'size': 8, 'weight': 'semibold'})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-----------------------------------------------------
#prep
l_update=l_update+['Unknown']
l_keep=[i for i in ada0.obs['predict'].astype('str').unique() if i not in l_update]

#update markers
for cell in l_keep: d_mark[cell]=d_mark0[cell]
with open(f'{fd_out}/mark.json', 'w') as f: json.dump(d_mark, f)  

#update ada
df=ada.obs.loc[:, ['prelabel', 'predict']].copy()
df=df.astype('str')
ada0.obs.loc[df.index, ['prelabel']]=df['prelabel']
ada0.obs.loc[df.index, ['predict']]=df['predict']
ada0.write(f'{fd_out}/data.h5ad')

#plot
f_out=f'{fd_out}/predict.png'
plt_dr(ada0, f_out, title='Predicted Cells')

