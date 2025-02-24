from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import json

#-----------------------------------------------------
l_sample=['Het', 'NonF', 'Fluc']
d_sample={'Het': 'Het NF', 'NonF': 'DE17.5 NF', 'Fluc': 'DE17.5 F'}

d_rep={
'raw': 'Raw',
'qc_cnt': 'Filter by Cell Counts',
'qc_gene': 'Filter by Gene Counts',
'qc_mt': 'Filter by Mt%',
'qc_dd': 'Dedoublet'
}

fd_out='./out/d00_plot_06_filter-cnt'
fd_in='./out/a02_clean_00_filter'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#-----------------------------------------------------
def get_df(f0, f1):
    df0=pd.read_csv(f'{fd_in}/{f0}.csv', index_col=0)
    df1=pd.read_csv(f'{fd_in}/{f1}.csv', index_col=0)
    #merge
    df0=df0.set_index('filter')
    df1=df1.set_index('filter')
    df=df0.merge(df1, left_index=True, right_index=True)
    #clean
    df['cnt']=df.sum(axis=1)
    df=df.reindex(['cnt'], axis=1)
    df['sample']=f0[:-1]
    df=df.reset_index()
    #rename
    df['filter']=df['filter'].replace(d_rep)
    return df


def plt_qc(df, f_out, title='Cell Counts After Sequential Filtering', sz=(8,3.5), cmap=['#666666', '#828282', '#999999', '#b0b0b0', '#ff0000']):
    #prep
    if sz is None: sz=(9, 0.1*df.shape[0])
    #plot
    fig, ax=plt.subplots(figsize=sz)
    ax=sns.barplot(x='sample', y='cnt', data=df, hue='filter', palette=cmap)
    #adjust
    sns.despine()
    ax.set_title(title, x=0.6, fontsize=22, weight='semibold', pad=10)
    plt.xlabel('')
    plt.ylabel('Counts', fontsize=15, labelpad=12, weight='semibold')
    plt.xticks(fontsize=14, weight='semibold')
    plt.legend(loc=(1.01, 0), frameon=False, prop={'weight': 'semibold', 'size':10})
    #save
    plt.tight_layout()
    plt.savefig(f_out, dpi=300)
    plt.close()
    return

#-----------------------------------------------------
#prep
df0=get_df('Het0', 'Het1')
df1=get_df('NonF0', 'NonF1')
df2=get_df('Fluc0', 'Fluc1')
df=pd.concat([df0, df1, df2])
df['sample']=df['sample'].replace(d_sample)

df.to_csv(f'{fd_out}/cnt.csv', index=False)

#plot
f_out=f'{fd_out}/cnt.png'
plt_qc(df, f_out)


