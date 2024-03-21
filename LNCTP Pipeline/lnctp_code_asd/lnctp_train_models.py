import pandas as pd
import numpy as np
import os
import tensorflow as tf
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.layers import Dense, BatchNormalization


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 21:05:49 2022
@ LNCTP testing models code
@ authors: Andy Yang, Jonathan Warrell
"""

class denseModel(tf.keras.Model):
    
    #initialize with customizable cell type, genes per cell type, dense nodes, activation function, and whether to use norm
    def __init__(self, n_param, dense_nodes, activation,batch_norm = True):
        
        super(denseModel,self).__init__()
        
        self.batch_norm = batch_norm
        
        #buildingh layers
        
        self.dense1 = Dense(dense_nodes,input_shape = (n_param,), activation = activation)
        
        if batch_norm:
            self.norm1 = BatchNormalization()
            
        self.dense2 = Dense(dense_nodes, activation = activation)
        
        if batch_norm:
            self.norm2 = BatchNormalization()
            
        self.out = Dense(1, activation = 'sigmoid')
        
    def call(self,x):
        #feed forward layers
        x = self.dense1(x)
        if self.batch_norm:
            x = self.norm1(x)
        x = self.dense2(x)
        if self.batch_norm:
            x = self.norm2(x)
        
        return self.out(x)
    
class sparseModel(tf.keras.Model):
    #same as other model but with embedding layer nodes for customization
    def __init__(self, embed_nodes, dense_nodes, activation, cell_types = 8, params =33, batch_norm = True):
        
        super(sparseModel,self).__init__()
        
        self.batch_norm = batch_norm
        
        #build layers with mask
        self.embed = Dense((cell_types+1)*embed_nodes,input_shape = (cell_types*params,), activation = activation,kernel_constraint = custom_mask(embed_nodes = embed_nodes))
        
        if batch_norm:
            self.norm1 = BatchNormalization()
            
        self.dense1 = Dense( dense_nodes, activation = activation)
        
        if batch_norm:
            self.norm2 = BatchNormalization()
            
        self.dense2 = Dense(dense_nodes, activation = activation)
        
        if batch_norm:
            self.norm3 = BatchNormalization()
            
        self.out = Dense(1, activation = 'sigmoid')
        
    def call(self,x):
        #feed forward
        x = self.embed(x)
        if self.batch_norm:
            x= self.norm1(x)
        x = self.dense1(x)
        if self.batch_norm:
            x = self.norm2(x)
        x = self.dense2(x)
        if self.batch_norm:
            x = self.norm3(x)
        
        return self.out(x)

#custom sparse layer mask with cell types, genes per cell type, and embedding node size
class custom_mask(tf.keras.constraints.Constraint):
    def __init__(self, embed_nodes, cell_types = 8, params =33):
        self.cell_types = cell_types
        self.params = params
        self.embed_nodes = embed_nodes
    
    def __call__(self, W):
        indices = []
        values = []
    
        #build sparse structure
        for i in range(self.cell_types):
            for j in range(self.params):
                for k in range(self.embed_nodes):
                    indices.append([i*self.params + j,self.embed_nodes*i + k])
                    values.append(1)
        i = self.cell_types
        for j in range(49):
            for k in range(self.embed_nodes):
                indices.append([i*self.params + j, self.embed_nodes*i + k])
                values.append(1)
    
        #create sparse structure
        sparse = tf.sparse.SparseTensor(indices = indices, values = values, dense_shape = [self.cell_types*self.params+49 , (self.cell_types+1)*self.embed_nodes])
    
        return tf.math.multiply(tf.cast(tf.sparse.to_dense(sparse),tf.float32),W)
    
    
        
        
class shallowSparseModel(tf.keras.Model):
    #same as other model but with embedding layer nodes for customization
    def __init__(self, embed_nodes, activation, cell_types = 8, params =33,  batch_norm = True):
        
        super(shallowSparseModel,self).__init__()
        
        self.batch_norm = batch_norm
        
        #build layers with mask
        self.embed = Dense((cell_types+1)*embed_nodes,input_shape = (cell_types*params,), activation = activation,kernel_constraint = custom_mask(embed_nodes = embed_nodes))
        
        if batch_norm:
            self.norm1 = BatchNormalization()
            
        self.out = Dense(1, activation = 'sigmoid')
        
    def call(self,x):
        #feed forward
        x = self.embed(x)
        if self.batch_norm:
            x= self.norm1(x)
        
        return self.out(x)


models = ['Dense', 'Deep Sparse','Shallow Sparse']
folds = [1,2,3,4,5,6,7,8,9,10]
dense_nodes = [5,10,20]
embed_nodes = [3,5,10]
learning_rates = [0.01,0.001,0.0001]
batch_sizes = [25,50,100,250]

accuracy_out = []

for fold in folds:

    X_train = pd.read_csv('data/train_'+str(fold)+'.csv',header = 0,index_col = 0)
    X_val = pd.read_csv('data/val_'+str(fold)+'.csv',header = 0,index_col = 0)
    X_test = pd.read_csv('data/test_'+str(fold)+'.csv', header =  0,index_col = 0)
    Y_train = pd.read_csv('data/train_trait_'+str(fold),index_col = 0,header = None)
    Y_val = pd.read_csv('data/val_trait_'+str(fold),index_col = 0,header = None)
    Y_test = pd.read_csv('data/test_trait_'+str(fold),index_col = 0,header = None)


    X_train = X_train.to_numpy()
    X_test = X_test.to_numpy()
    X_val = X_val.to_numpy()

    fold_loss = []
    spec_models = []
    for mod in models:
        for embed_node in embed_nodes:
            for dense_node in dense_nodes:
                for lr in learning_rates:
                    for batch_size in batch_sizes:

                        if mod=="Dense":
                            model = denseModel(n_param = 321, dense_nodes = dense_node, activation = 'linear')


                        elif mod=="Deep Sparse":
                            model = sparseModel(embed_nodes = embed_node, dense_nodes = dense_node, activation = 'linear')


                        elif mod=='Shallow Sparse':
                            model = shallowSparseModel(embed_nodes = embed_node, activation = 'linear')
                    
                        model.compile(loss='binary_crossentropy', optimizer=tf.keras.optimizers.Adam(learning_rate=lr), metrics=['accuracy'])
                
                        early_stopping = EarlyStopping(monitor='val_accuracy',patience=8, min_delta=0.001, mode='max')


                        model.fit(X_train,Y_train,epochs = 100, batch_size = batch_size, validation_data=(X_val,Y_val),callbacks=[early_stopping])

        
                        results = model.evaluate(X_test, Y_test, batch_size = 128)
                        accuracy = results[1]     
                        fold_loss.append(accuracy)
                        
                        if mod=="Dense":
                            spec = 'Dense_denseNodes'+str(dense_node)+'_embedNodes+'+str(embed_node)+'_lr'+str(lr)+'_batchsize'+str(batch_size)
                            spec_models.append('Dense_denseNodes'+str(dense_node)+'_embedNodes+'+str(embed_node)+'_lr'+str(lr)+'_batchsize'+str(batch_size))


                        elif mod=="Deep Sparse":
                            spec = 'DeepSparse_denseNodes'+str(dense_node)+'_embedNodes+'+str(embed_node)+'_lr'+str(lr)+'_batchsize'+str(batch_size)
                            spec_models.append('DeepSparse_denseNodes'+str(dense_node)+'_embedNodes+'+str(embed_node)+'_lr'+str(lr)+'_batchsize'+str(batch_size))


                        elif mod=='Shallow Sparse':
                            spec = 'ShallowSparse_denseNodes'+str(dense_node)+'_embedNodes+'+str(embed_node)+'_lr'+str(lr)+'_batchsize'+str(batch_size)
                            spec_models.append('ShallowSparse_denseNodes'+str(dense_node)+'_embedNodes+'+str(embed_node)+'_lr'+str(lr)+'_batchsize'+str(batch_size))

                        #preds = model.predict(X_test)
                        #pred_df = pd.DataFrame(preds)
                        #pred_df.to_csv('preds/' + spec + '_' + str(fold) + '.csv')

                        model.save('models/' + spec + '_' + str(fold))

    losses = pd.Series(data = fold_loss, index = spec_models)
    accuracy_out.append(losses)
    accuracy_out2 = pd.concat(accuracy_out,axis = 1)
    accuracy_out2 = accuracy_out2.round(3)
    accuracy_out2.to_csv('model_accuracy_summary.csv')

