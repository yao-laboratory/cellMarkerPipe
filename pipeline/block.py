#!/usr/bin/env python
# coding: utf-8

#import pipeline
#import pipeline.os as os
#import pipeline.time as time
#import pipeline.pd as pd
from pipeline import *
#print(os)
code_dir = PATH_TO_SOURCE + "/" + "src/"


# step1: preprocess
def preprocess(work_dir, data_dir, nvariable=2000, max_RNA = 2500, min_RNA = 5, max_mt = 5, nPCA=10, resolution=0.5, algorithm=1, Cluster=False, know_marker=False, keep_known=False):

    # preprocess
    command = (
        "Rscript "
        + code_dir
        + "preprocess.r "
        + work_dir
        + " "
        + data_dir
        + " "
        + str(nvariable)
        + " "
        + str(max_RNA)
        + " "
        + str(min_RNA)
        + " "
        + str(max_mt)
        + " "
        + str(nPCA)
        + " "
        + str(resolution)
        + " "
        + str(algorithm)
        + " "
        + str(Cluster)
        + " "
        + str(know_marker)
        + " "
        + str(keep_known)
        + " > stat_preprocess"
    )

    os.system(command)


    return



# step 2: Selection


def selection(work_dir, data_dir="", method="de", n_marker=10,  **kwarg):

    # map the inputs to the function blocks
    options = {
        "de": diff_express,  # by group
        "scG": scGenefit,  # not by group
        "Com": Comet,  # by group, by neighbouring number
        "SC": SCmarker,  # not by group
        "cos": COSGmarker,  # by group
        "fst": FEASTmarker,  # not by group
        "sc3": SC3_diff, # by group
        "hv": high_variable # by group
    }

    options[method](work_dir=work_dir, data_dir=data_dir, n_marker=n_marker,  **kwarg)
    return



# step 3: Evaluation

# supervised evaluation
def evaluation(work_dir, nPCA=10,  resolution=0.5, algorithm=1, know_marker=False, data_dir=""):

    # cluster again
    #command = "Rscript " + code_dir + "re-cluster.r " + work_dir  + " " + str(nPCA) + " " + str(resolution) + " " + str(algorithm)
   
    #os.system(command)

    # compare y_predict with y_true
    command = "python " + code_dir + "evaluation.py " + work_dir + " " + str(know_marker) + " " + data_dir + " > stat_evaluation"

    
    os.system(command)

    return



# seurat differential analyses method
def diff_express(work_dir, data_dir="", n_marker=10):
    # n is the number of marker for each cluster

    command = (
        "Rscript "
        + code_dir
        + "diff_express.r "
        + work_dir
        + " "
        + str(n_marker)
        + " > stat_selection"
    )
   
    os.system(command)

    return



# SC3 method
# use envrionment yingluR4.1
def SC3_diff(work_dir,data_dir="", n_marker=10):
    # n is the number of marker for each cluster

    command = (
        "Rscript "
        + code_dir
        + "diff_SC3.r "
        + work_dir
        + " "
        + str(n_marker)
        + " > stat_selection"
    )
    os.system(command)

    return



# scGenefit method
def scGenefit(
    work_dir, data_dir="",input_format="10X", n_marker=10, method="centers", epsilon= 1, redundancy=0.25):
    # code directory

    # input files
    data_path = work_dir + "/data/norm_data_high_variable.csv"
    label_path = work_dir + "/data/cluster_labels.csv"
    name_path = work_dir + "/cluster/tabnames.csv"

    command = (
        "python "
        + code_dir
        + "/scGfit.py "
        + work_dir
        + " "
        + input_format
        + " "
        + str(n_marker)
        + " "
        + str(epsilon)
        + " > stat_selection"
    )

    os.system(command)

    return


# Comet method
# use envrionment py36
def Comet(work_dir, data_dir="", n_marker=10, vis=False, if_pair=1, others=""):
    # data_dir is the directory that users have their own 3 files "tabcluster.txt  tabmarker.txt  tabvis.txt"
    # otherwise, the user needs cluster first, and the input files will be prduced by the last "preprocess" step
    # vis: wehter visulize with Comet, if vis==True, then need to provide a tabvis.txt file under data_dir, or work_dir/cluster
    # if_pair: choose the markers from the paired result(1); or single result (0)
    # n_marker: the number of markers (or pairs) to choose
    # others: other parameters for Comet, should be a string with the format described by Comet manual
    # this job needs python 36. Should be run under environment Py36
    # code directory
    if data_dir == "":
        # tried to use normalized data, error pop up
        tabmarker_file = work_dir + "/data/scale_data_high_variable.csv"
        tabcluster_file = work_dir + "/data/cluster_labels.csv"
        tabvis_file = work_dir + "/data/tabvis.txt"

        command = (
            "Comet "
            + tabmarker_file
            + " "
            + tabvis_file
            + " "
            + tabcluster_file
            + " marker/ "
            + others
        )
        if not vis:
            with open(tabvis_file, "w") as fp:
                pass
            command = (
                "Comet "
                + tabmarker_file
                + " "
                + tabvis_file
                + " "
                + tabcluster_file
                + " marker/ -skipvis True "
                + others
            )

    else:
        tabmarker_file = data_dir + "/scale_data_high_variable.csv"
        tabcluster_file = data_dir + "/cluster_labels.csv"
        tabvis_file = data_dir + "/tabvis.txt"
        command = (
            "Comet "
            + tabmarker_file
            + " "
            + tabvis_file
            + " "
            + tabcluster_file
            + " marker/ "
            + others
        )
        if not vis:
            with open(tabvis_file, "w") as fp:
                pass
            command = (
                "Comet "
                + tabmarker_file
                + " "
                + tabvis_file
                + " "
                + tabcluster_file
                + " marker/ -skipvis True "
                + others
            )

    os.system(command)
    # get top # of markers
     # compare y_predict with y_true
    command = ("python " 
               + code_dir 
               + "/comet_result.py " 
               + work_dir
               + " "
               + str(n_marker)
               + " "
               + str(if_pair)
               + " > stat_comet_result")

    os.system(command)
    
    return



# SCmarker method
def SCmarker(work_dir, data_dir="", k=100, n=10, n_marker=10):

    command = (
        "Rscript "
        + code_dir
        + "scmarker.r "
        + work_dir
        + " "
        + data_dir
        + " "
        + str(k)
        + " "
        + str(n)
        + " "
        + str(n_marker)
        + " > stat_selection"
    )

    os.system(command)

    return


# COSG method
def COSGmarker(work_dir, data_dir="",n_marker=10, mu=1):
    # n is the number of marker for each cluster

    command = (
        "Rscript "
        + code_dir
        + "cosg_marker.r "
        + work_dir
        + " "
        + str(n_marker)
        + " "
        + str(mu)
        + " > stat_selection"
    )
   
    os.system(command)

    return

# FEAST method
def FEASTmarker(work_dir,data_dir="", n_marker=10):
    # n is the number of marker for each cluster

    command = (
        "Rscript "
        + code_dir
        + "FEAST_marker.r "
        + work_dir
        + " "
        + str(n_marker)
        + " > stat_selection"
    )

    os.system(command)

    return


# high vaiable genes method
def high_variable(work_dir, data_dir="", n_marker=10):

    command = (
        "Rscript "
        + code_dir
        + "high_variable.r "
        + work_dir
        + " "
        + data_dir
        + " "
        + str(n_marker)
        + " > stat_selection"
    )

    os.system(command)

    return
