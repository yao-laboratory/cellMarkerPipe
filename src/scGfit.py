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

n_marker = sys.argv[2]
n_marker = int(n_marker)
print("The number of markers is:")
print(n_marker)

input_format = sys.argv[3]
print("The input format is:")
print(input_format)

epsilon = sys.argv[4]
epsilon = float(epsilon)
print("Epsilon is chosen as:")
print(epsilon)
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

# optimize parameter epsilons
#create data
#optional: generate new data to evaluate performance

#X_train, y_train=create_data(f_groups, cell_types, 0.5, n, [10,10,100])
#X_test, y_test=create_data(f_groups, cell_types, 0.5, n, [10,10,100])

#res=optimize_epsilon(X_train, y_train, X_test, y_test, m, method='centers', fixed_parameters={'redundancy': 0, 'max_constraints': 100},bounds=[(0.2 , 10)], x0=[1], max_fun_evaluations=25, n_experiments=10,clf=clf, verbose=True)

#print("best epsilon:", res[0], "accuracy:", res[1])
# find markers
start = time.time()
markers= get_markers(data, labels, n_marker*len(Y_t), method=centers, epsilon = epsilon,redundancy=0.25)
end = time.time()
print("Excuting time:%s second"%(end-start))
print(genes[markers])
all_markers = pd.DataFrame (genes[markers], columns = ['markers']) 
select_dir = work_dir + "/marker"

if not os.path.isdir(select_dir):
    os.mkdir(select_dir)

all_markers_dir =  select_dir + "/all_genes.txt" 
all_markers.to_csv(all_markers_dir, header = True, sep='\t')
# to get the marker gene for each cluster, the average expression of each selected genes for each cluster is going to be calculated
# then the top ones is going to be assignend to the cluster
# subset the dataset with the selected markers
data= data[:,markers]
# res save the mean expression of each gene in each cluster. dimension clusters * seleceted genes
res = pd.DataFrame(columns=genes[markers], index=Y_t)
print(res)
print(res.shape)

for clust in Y_t:
    res.loc[clust] = data[Y_true==clust,:].mean(0)
# marker_per group save the markers of each group. dimension: clusters * marker genes
#marker_per_group = pd.DataFrame( index= Y_t, columns=range(n_marker))
marker_per_group = pd.DataFrame( index= Y_t, columns=["Marker"])
marker_per_group.index.name = "Cluster"
# unique_markers save all unique markers of all clusters
unique_markers = pd.Index([])
for clust in Y_t:
    markers = res.loc[clust].sort_values(ascending = False)[:n_marker,].index
    #marker_per_group.loc[clust] = markers
    marker_per_group.loc[clust] = ", ".join(markers)
    unique_markers = unique_markers.append(markers)

unique_markers = unique_markers.unique()
df_markers = pd.DataFrame (unique_markers, columns = ['markers'])

# save the result
result_dir = select_dir + "/marker_genes.txt"
result_per_group_dir = select_dir + "/marker_gene_per_group.csv"
df_markers.to_csv(result_dir, header = True, sep='\t')
#marker_per_group.to_csv(result_per_group_dir, header = True, sep='\t')
marker_per_group.to_csv(result_per_group_dir, header = True, sep=',')

