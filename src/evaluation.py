#!/usr/bin/env python
# coding: utf-8

from tools import *

import sys
import csv
import pandas as pd
import os
import re
# total arguments
n = len(sys.argv)
#print("Total arguments passed:", n)
 
# Arguments passed
#print("\nName of Python script:", sys.argv[0])
 
#print("\nArguments passed:", end = " ")
#for i in range(1, n):
#    print(sys.argv[i], end = " ")
     
# Addition of numbers
#Sum = 0
# Using argparse module
#for i in range(1, n):
#    Sum += int(sys.argv[i])
     
#print("\n\nResult:", Sum)

work_dir = sys.argv[1]
print("Work directory is:")
print(work_dir)
know_marker = sys.argv[2]
print("Do you have a file stored known marker genes:")
print(know_marker)
if know_marker == "True":
    data_dir = sys.argv[3]
    print("Where the Known_marker.csv saved?")
    print(data_dir)

cluster_dir = work_dir + "/data/cluster_labels.csv"
re_cluster_dir = work_dir + "/re-cluster/re_cluster_labels.csv"

print("Need to find the cluster information here:")
print(cluster_dir)
print("To calculate evaluation index, the cells are re-clustered here:")
print(re_cluster_dir)
#labels_true = [0, 0, 1, 1, 1, 1]
#labels_pred = [0, 0, 2, 2, 3, 3]
# open the file in the write mode
eval_dir = work_dir + "/evaluation"
if not os.path.isdir(eval_dir):
    os.mkdir(eval_dir)

#labels_true = pd.read_csv(cluster_dir)
labels_true = pd.read_csv(cluster_dir, sep = '\t', header=None)
labels_pred = pd.read_csv(re_cluster_dir)

print(labels_true.iloc[:,1].values)
print(labels_pred.iloc[:,1].values)
# matrix plot of the marker genes

Y_true = labels_true.iloc[:,1].values
Y_pred = labels_pred.iloc[:,1].values
# calculate the evaluation values
Y_t, indices = np.unique(Y_true, return_inverse=True)

eval_result = eval_cluster(indices, Y_pred)
eval_result
print(eval_result)

if know_marker == "True":
    Known_marker_dir = data_dir + "/Known_marker.csv"
    # import markers genes of well known
    colnames = ["Cluster", "Markers"]
    Know_marker_genes = pd.read_csv(Known_marker_dir, header= None, sep=",", names=colnames)
    dfk = Know_marker_genes.transpose()
    cdfk = dfk.rename(columns=dfk.loc["Cluster"].astype(str))
    for index in cdfk.columns:
        cdfk[index]["Markers"]=re.split(', ',cdfk[index]["Markers"])
        cdfk[index]["Markers"]=[i.upper() for i in cdfk[index]["Markers"]]
    cdfk=cdfk.drop(labels="Cluster", axis=0)
    print(cdfk)
    
    # import marker genes found by the method
    markerfile= work_dir + "/marker/marker_gene_per_group.csv"
    marker_genes = pd.read_csv(markerfile)
    print("Found the markers")
    top_rows = marker_genes.head()
    print(top_rows)

    df = marker_genes.transpose()
   # print(df)
    cdf = df.rename(columns=df.loc["Cluster"].astype(str))
   # print(cdf)
    cdf=cdf.drop(labels="Cluster", axis=0)
    print("reshape the marker dataframe")
    top_rows = cdf.head()
    print(top_rows)
#    cdf=cdf.sort_index(axis=1)
    for index in cdf.columns:
        print(cdf[index]["Marker"])
        print(pd.__version__)
        cdf[index]["Marker"]=cdf[index]["Marker"].split(", ")
        cdf[index]["Marker"]=[i.upper() for i in cdf[index]["Marker"]]
    cdf = cdf.reset_index(drop = True)
    marker_genes_dict = cdf.to_dict('records')[0]
    df_found = pd.DataFrame(columns=cdfk.columns, index=["found_markers", "precision_per_group","recal_per_group", "tot_precision", "tot_recal"])
    tot_found = 0
    tot_len_b = 0
    tot_len_a = 0
    for index in cdfk.columns:
        a = cdfk[index][0]
        b = cdf[index][0]
        found_markers = set(a).intersection(set(b))
        df_found[index]["found_markers"]=found_markers
        df_found[index]["precision_per_group"] = len(found_markers)/len(b)
        df_found[index]["recal_per_group"] = len(found_markers)/len(a)
        tot_found = tot_found + len(found_markers)
        tot_len_b = tot_len_b + len(b)
        tot_len_a = tot_len_a + len(a)
    df_found.iloc[3,0] = tot_found/tot_len_b
    df_found.iloc[4,0] = tot_found/tot_len_a
    save_found_dir = eval_dir + "/precision_recall.csv"
    df_found.to_csv(save_found_dir, header = True, sep='\t')

# save the result to file
result_dir = eval_dir + "/result.csv"
eval_result.to_csv(result_dir, header = True)
#f = open(result_dir, 'w')

# create the csv writer
#writer = csv.writer(f)

# write a row to the csv file
#writer.writerow(eval_result)

# close the file
#f.close()
