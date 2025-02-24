from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import tensorflow.keras as keras
import random
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
import seaborn as sns

#-----------------------------------------------------
k=3
n_gene=1000

fd_out='./out/a01_dd_01_r0'
fd_in='./out/a01_dd_00_prep'

######################################################
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_fname=list(Path(fd_in).glob('*.h5ad'))

#-----------------------------------------------------
class seq_dense: 
    def __init__(self, n_layer=1, lr=0.0001, r=0.1):
        self.layer=n_layer
        self.lr=lr
        self.r=r

    def build(self, n_in, n_out):
        #model
        reg=keras.regularizers.L1L2(l1=0.1, l2=0.1)
        model=keras.models.Sequential()
        model.add(keras.layers.Input(n_in,))
        for i in range(self.layer):
            model.add(keras.layers.Dropout(self.r))
            model.add(keras.layers.Dense(n_in*2, activation='swish', kernel_initializer='he_normal', kernel_regularizer=reg))
        model.add(keras.layers.Dense(n_out, activation='softmax'))
        print(model.summary())
        #compile
        opt=keras.optimizers.legacy.Adam(learning_rate=self.lr)
        model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['categorical_accuracy'])
        return model

def sim_dd(X):
    #prep
    ar_l=np.sum(X, axis=1)  #all library size
    n=X.shape[0] 
    #generate multiplets copies (total n*n_sim simulate cells)
    ar_mul=np.random.poisson(lam=2, size=n)
    ar_mul=np.where(ar_mul<=2, 2, ar_mul)  #at least doubles
    #sim
    l_sim=[]
    for i in ar_mul:
        idx=np.random.randint(n, size=i)
        ar_sim=X[idx, :].copy()
        n_lib=np.sum(ar_sim, axis=1).max()  #largest library size in the picked row, for normalization
        ar_sim=ar_sim.sum(axis=0)
        #norm
        n_lib=np.random.choice(ar_l[ar_l>=n_lib])
        sum_sim=ar_sim.sum()
        ar_sim=ar_sim*n_lib/sum_sim
        ar_sim=ar_sim.astype('int')
        l_sim.append(ar_sim)
    X_sim=np.stack(l_sim)
    return X_sim


#-----------------------------------------------------
for fname in l_fname:
    #prep
    sample=Path(fname).stem
    ada=sc.read(fname)
    ada_save=ada.copy()
    df=ada.to_df()
    
    #non zero genes
    dfi=df.T.copy()
    dfi['pct']=(dfi>0).mean(axis=1)
    dfi=dfi.sort_values('pct', ascending=False)
    l_gene=dfi.index.tolist()[:n_gene]

    #split train
    l_df=np.array_split(df, k)
    l_pred=[]
    for i, dfi in enumerate(l_df):
        #simulate cells (sim on all genes)
        df_sim=pd.DataFrame(sim_dd(dfi.values), columns=dfi.columns)
        df_sim=df_sim.reindex(l_gene, axis=1)
        X_sim=df_sim.values
        #original cells
        df_ori=dfi.reindex(l_gene, axis=1)
        X_ori=df_ori.values
        #train X
        X_train=np.log(np.vstack([X_ori, X_sim])+1)
        #train y
        y0=np.zeros(X_ori.shape[0])  #original is False (not doublets)
        y1=np.ones(X_sim.shape[0]) 
        y_train=np.concatenate([y0, y1])
        X, y=shuffle(X_train, y_train, random_state=42)
        y=keras.utils.to_categorical(y)
        #test
        idx=[j for j in df.index if j not in dfi.index]
        df_test=df.reindex(idx).reindex(l_gene, axis=1)
        X_test=np.log(df_test.values+1)
        #model
        model=seq_dense()        
        mod=model.build(n_in=X.shape[1], n_out=2)
        mod.fit(X, y, epochs=10, batch_size=32, validation_split=0.1)#, callbacks=[cb_checkpoint])
        #predict (doublet prob)
        df_pred=pd.DataFrame(mod.predict(X_test)[:, 1], index=df_test.index, columns=[f'mod_{i}'])
        l_pred.append(df_pred)
        
    #concat pred (score is doublet prob)
    df=pd.concat(l_pred, join='outer', axis=1)
    df['dd_score']=df.mean(axis=1)
    df=df.reindex(['dd_score'], axis=1)

    #add to ada
    ada.obs=ada.obs.merge(df, left_index=True, right_index=True)
    ada.write(f'{fd_out}/{sample}.h5ad')

    print(ada.obs)
