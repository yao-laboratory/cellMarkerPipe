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
unique_markers = np.array([],dtype=object)
marker_per_group = pd.DataFrame(index=Y_t, columns=["Marker"])
marker_per_group.index.name = "Cluster"
for group in Y_t:

    if (if_pair):
        marker_group_file = output_dir + "/cluster_" + str(group) + "_pair_final_ranking.csv"
        marker_file = pd.read_csv(marker_group_file, header=0)
        markers = []
        for i in marker_file.index:
            if(len(markers) >= n_marker):
                break
            marker1 = marker_file["gene_1"][i].split("_")[0]
            if (marker1 not in markers):
                markers+= [marker1]
            if(len(markers) >= n_marker):
                break
            marker2 = marker_file["gene_2"][i].split("_")[0]
            if (marker2 not in markers):
                markers+= [marker2]
        marker_per_group.loc[group] = ", ".join(markers)
        unique_markers = np.unique(np.concatenate((unique_markers, markers)))


    else:
        marker_group_file = output_dir + "/cluster_" + str(group) + "_singleton_all_ranked.csv"
        marker_file = pd.read_csv(marker_group_file, header=0)
        markers = marker_file.iloc[:,0].str.split(pat="_",expand=True)[0].unique()[0:n_markers]
        unique_markers = np.unique(np.concatenate((unique_markers, markers)))
        marker_per_group.loc[group] = ", ".join(markers)


# save the result
df_markers = pd.DataFrame (unique_markers, columns = ['markers'])
#result_dir = select_dir + "/marker_genes.txt"
#df_markers.to_csv(result_dir, header = True, sep='\t')
#result_dir = select_dir + "/marker_gene_per_group.csv"
#df_marker_group.to_csv(result_dir, header = True, sep='\t')

result_dir = select_dir + "/marker_genes.txt"
result_per_group_dir = select_dir + "/marker_gene_per_group.csv"
print(result_per_group_dir)
df_markers.to_csv(result_dir, header = True, sep='\t')
marker_per_group.to_csv(result_per_group_dir, header = True, sep=',')
