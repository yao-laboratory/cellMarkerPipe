#!/usr/bin/env python
# coding: utf-8

# ### Contingency Matrix 

# In[4]:


from sklearn.metrics.cluster import contingency_matrix 
from sklearn.metrics.cluster import adjusted_rand_score 
from sklearn.metrics.cluster import normalized_mutual_info_score 
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics import silhouette_score
from sklearn.metrics import pairwise_distances
from sklearn.metrics import jaccard_score
from sklearn.metrics import f1_score
from sklearn import datasets
import pandas as pd
import numpy as np

def purity_score(y_true, y_pred): 
    # compute contingency matrix (also called confusion matrix) 
    c_matrix = contingency_matrix(y_true, y_pred) 
    # return purity 
    return np.sum(np.amax(c_matrix, axis=0)) / np.sum(c_matrix)  

def jaccard_similarity(list1, list2): 
    s1 = set(list1) 
    s2 = set(list2) 
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))


def eval_cluster(labels_true, labels_pred):
    ari = adjusted_rand_score(labels_true, labels_pred)
    jaccard = jaccard_score(labels_true, labels_pred, average='macro') 
    purity = purity_score(labels_true, labels_pred)
    nmi = normalized_mutual_info_score(labels_true, labels_pred)
    fmi = fowlkes_mallows_score(labels_true, labels_pred)
    f1_mac = f1_score(labels_true, labels_pred, average='macro')
    f1_mic = f1_score(labels_true, labels_pred, average='micro')

    column_names = ['ari', 'jaccard', 'purity', 'nmi', 'fmi', 'f1_macro', 'f1_micro']
    matrix = np.reshape([ari, jaccard, purity, nmi, fmi, f1_mac, f1_mic], (1,7))
    
    res = pd.DataFrame(matrix, columns=column_names)
    return res

