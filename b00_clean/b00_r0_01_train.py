from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import json
import tensorflow.keras as keras
import matplotlib.pyplot as plt
import seaborn as sns

#-----------------------------------------------------
#n_input=1000
min_p=0.5

fd_out='./out/b00_r0_01_train'
fd_in='./out/b00_r0_00_mark'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)

#gene
dfg=pd.read_csv(f'{fd_in}/gm_filtered.csv', index_col=0)
l_gene=dfg.index.tolist()

#ada
ada=sc.read(f'{fd_in}/data.h5ad')
ada_save=ada.copy()
ada=ad.AnnData(ada.raw.X, obs=ada.obs, var=ada.raw.var)
ada.obs['prelabel']=ada.obs['prelabel'].astype('str')

#-----------------------------------------------------
class seq_dense:
    def __init__(self, n_layer=1, r=0.1, lr=0.0001, activation='swish', l1=0, l2=0):
        self.layer=n_layer
        #self.node=n_gene*2
        self.r=r
        self.lr=lr
        self.acti=activation
        self.l1=l1
        self.l2=l2
    def build(self, n_in, n_out):
        #model
        reg=keras.regularizers.L1L2(l1=self.l1, l2=self.l2)
        model=keras.models.Sequential()
        model.add(keras.layers.Input(n_in,))
        for i in range(self.layer):
            model.add(keras.layers.Dropout(self.r))
            model.add(keras.layers.Dense(n_in*2, activation=self.acti, kernel_initializer='he_normal', kernel_regularizer=reg))
        model.add(keras.layers.Dense(n_out, activation='softmax'))
        print(model.summary())
        #compile
        opt=keras.optimizers.legacy.Adam(learning_rate=self.lr)
        model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['categorical_accuracy'])
        return model


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
#input df (all cell) with filtered genes
df=ada.to_df().reindex(l_gene, axis=1)

#model
l_cell=[i for i in ada.obs['prelabel'].unique() if i!='TBD']
model=seq_dense()
model=model.build(df.shape[1], len(l_cell))

#Xy
df_train=ada.obs.query('prelabel in @l_cell')
X=df.reindex(df_train.index).values
y=df_train['prelabel'].replace(l_cell, range(len(l_cell))).astype('int').values.flatten()
y=keras.utils.to_categorical(y)

#train
model.fit(X, y, epochs=10, batch_size=32, validation_split=0.1, callbacks=[])

#pred
df_pred=pd.DataFrame(model.predict(df.values), columns=l_cell, index=df.index)
df_pred['predict']=df_pred.idxmax(axis=1)
df_pred.loc[df_pred.max(axis=1, numeric_only=True)<min_p, ['predict']]='Unknown'
df_pred.to_csv(f'{fd_out}/predict.csv')

#update ada
df_pred['prelabel']=ada_save.obs['prelabel'].astype('str')
df_pred['predict']=np.where(df_pred['prelabel']=='TBD', df_pred['predict'], df_pred['prelabel'])
df_pred=df_pred.loc[:, ['predict']].copy()
ada_save.obs=ada_save.obs.merge(df_pred, left_index=True, right_index=True)
ada_save.write(f'{fd_out}/data.h5ad')

#-----------------------------------------------------
#plot
ada=sc.read(f'{fd_out}/data.h5ad')
f_out=f'{fd_out}/pred.png'
plt_dr(ada, f_out, title='Predicted Cells')
