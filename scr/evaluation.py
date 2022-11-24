#!/usr/bin/env python
# coding: utf-8

from tools import *

import sys
import csv
import pandas as pd
import os
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

cluster_dir = work_dir + "/data/cluster_labels.csv"
re_cluster_dir = work_dir + "/re-cluster/re_cluster_labels.csv"

print("Need to find the cluster information here:")
print(cluster_dir)
print("To calculate evaluation index, the cells are re-clustered here:")
print(re_cluster_dir)
#labels_true = [0, 0, 1, 1, 1, 1]
#labels_pred = [0, 0, 2, 2, 3, 3]

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

# save the result to file

# open the file in the write mode
eval_dir = work_dir + "/evaluation"

if not os.path.isdir(eval_dir):
    os.mkdir(eval_dir)
result_dir = eval_dir + "/result.csv"
eval_result.to_csv(result_dir, header = True)
#f = open(result_dir, 'w')

# create the csv writer
#writer = csv.writer(f)

# write a row to the csv file
#writer.writerow(eval_result)

# close the file
#f.close()
