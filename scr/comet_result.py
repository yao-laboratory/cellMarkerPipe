#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import pandas as pd

# input variables
work_dir = sys.argv[1]
print("Work directory is:")
print(work_dir)


n_marker = sys.argv[2]
n_marker = int(n_marker)
print("The number of selected markers (pairs):")
print(n_marker)


if_pair = sys.argv[3]
if_pair = int(if_pair)
print("Select pair genes?")
if (if_pair):
    print("yes!")
    
else:
    print("no!")


# read the cluster file to loop each cluster
tabcluster_file = work_dir + "/data/cluster_labels.csv"
clusters = pd.read_csv(tabcluster_file, sep = '\t', header=None)
#print(clusters.shape)

output_dir = work_dir + "/marker/data"
clusters_np = clusters.iloc[:,1].values
Y_t, indices = np.unique(clusters_np, return_inverse=True)


select_dir = work_dir + "/marker"
# loop each cluster to read the result of COMET
marker_all = np.array([],dtype=object)
df_marker_group = pd.DataFrame(columns = ['Cluster', 'Markers'])
for group in Y_t:

    if (if_pair):
        marker_group_file = output_dir + "/cluster_" + str(group) + "_pair_final_ranking.csv"
        marker_file = pd.read_csv(marker_group_file, header=0)
        markers1 = marker_file.iloc[0:n_marker,1].str.split(pat="_",expand=True)[0].unique()
        markers2 = marker_file.iloc[0:n_marker,2].str.split(pat="_",expand=True)[0].unique()
        print(markers1)
        print(markers2)
        print(type(markers1))
        print(type(markers2))
        print(markers1.shape)
        print(markers2.shape)
        markers = np.concatenate((markers1, markers2),axis=0)
#        markers = np.concatenate(markers1, markers2)
        marker_all = np.unique(np.concatenate((marker_all, markers)))
       
        df_marker_group = df_marker_group.append({'Cluster':group, 'Markers':markers}, ignore_index=True)


    else:
        marker_group_file = output_dir + "/cluster_" + str(group) + "_singleton_all_ranked.csv"
        marker_file = pd.read_csv(marker_group_file, header=0)
        markers = marker_file.iloc[0:n_marker,0].str.split(pat="_",expand=True)[0].unique()
        marker_all = np.unique(np.concatenate((marker_all, markers)))

        df_marker_group = df_marker_group.append({'Cluster':group, 'Markers':markers})

# save the result
df_markers = pd.DataFrame (marker_all, columns = ['markers'])
result_dir = select_dir + "/marker_genes.txt"
df_markers.to_csv(result_dir, header = True, sep='\t')
result_dir = select_dir + "/marker_gene_per_group.csv"
df_marker_group.to_csv(result_dir, header = True, sep='\t')
