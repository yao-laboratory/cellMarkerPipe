#!/usr/bin/env python
# coding: utf-8

import sys
import os
import scanpy as sc
from scGeneFit.functions import *
import numpy as np
import time
np.random.seed(0)

from sklearn.neighbors import NearestCentroid
clf=NearestCentroid()
import pandas as pd

# input variables
work_dir = sys.argv[1]
print("Work directory is:")
print(work_dir)

input_format = sys.argv[2]
print("The input format is:")
print(input_format)

n_marker = sys.argv[3]
n_marker = int(n_marker)
print("The number of markers is:")
print(n_marker)

# read normalized matrix
if (input_format == "10X"):
# if it is 10x data, then read 10x data directly
    data_path=work_dir + "/data/10x"
    adata = sc.read_10x_mtx(
    data_path,  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=False)
    genes = adata.var_names
# scale the data
    sc.pp.scale(adata, max_value=10)
    data = adata.X
elif (input_format == "norm"):
# if the data format is matrix in csv, then read csv
# this section is saved only if users has a matrix which is already scaled, I did not provide scale here
# if the data is only log(x+1), then this algorithm is not able to identify enough genes
    data_path=work_dir + "/data/norm_data_high_variable.csv"
    data_df = pd.read_csv(data_path,sep='\t', index_col=0).rename_axis('cell',axis=1)
    data=data_df.to_numpy(dtype='float')
# need to add a scale algorithm here
    # this algorithm needs cell * gene
    data=np.transpose(data)
else:
# if the data format is scaled matrix, then read csv directly
# do not need to scale further
    data_path=work_dir + "/data/scale_data_high_variable.csv"
    data_df = pd.read_csv(data_path,sep='\t', index_col=0).rename_axis('cell',axis=1)
    data=data_df.to_numpy(dtype='float')
    # this algorithm needs cell * gene
    data=np.transpose(data)

# import cluster lable file
cluster_dir = work_dir + "/data/cluster_labels.csv"
labels_true = pd.read_csv(cluster_dir, sep = '\t', header=None)

# translate the string label into integer, then turn to dtype of unit8
Y_true = labels_true.iloc[:,1].values
Y_t, indices = np.unique(Y_true, return_inverse=True)
labels = indices.astype('uint8')

centers = "centers"
print(data.shape)
print(data)
print(labels.shape)
print(labels.dtype)

# find markers
start = time.time()
markers= get_markers(data, labels, n_marker, method=centers, redundancy=0.25)
end = time.time()
print("Excuting time:%s second"%(end-start))
print(genes[markers])

df_markers = pd.DataFrame (genes[markers], columns = ['markers'])

# save the result
select_dir = work_dir + "/marker"
if not os.path.isdir(select_dir):
    os.mkdir(select_dir)
result_dir = select_dir + "/marker_genes.txt"
df_markers.to_csv(result_dir, header = True, sep='\t')

