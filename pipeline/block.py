#!/usr/bin/env python
# coding: utf-8

#import pipeline
#import pipeline.os as os
#import pipeline.time as time
#import pipeline.pd as pd
from pipeline import *
#print(os)
code_dir = PATH_TO_SOURCE + "/" + "src/"
# function to check the status of the job
def stats(dep=""):
    command = "samtools flagstat output.sam > stats.txt"
    job_id = sbatch("stats", command, dep=dep)
    return job_id



# function to submitting jobs on hcc
def sbatch(
    job_name, command, work_dir= "", time=48, mem=60, tasks=1, partition="yaolab,batch", environment="gene_select", dep=""
):
    # jobs depend on the results of former jobs
    print("dep is:")
    print(dep)
    if dep != "":
        dep = "--dependency=afterok:{} --kill-on-invalid-dep=yes ".format(dep)
        # print(dep)

    # prepare a submitting file
    print("The current directory is:")
    print(os.getcwd())
    print("The work directory is:")
    print(work_dir)
    file = open(work_dir + "/sub.sh", "w")
    headline = "#!/bin/bash\n#SBATCH --time={}:00:00\n#SBATCH --mem-per-cpu={}G\n#SBATCH --job-name={}\n#SBATCH --error={}.err\n#SBATCH --output={}.out\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node={}\n#SBATCH --partition={}\n#SBATCH {}\n".format(
        time, mem, job_name, job_name, job_name, tasks, partition, dep
    )
    ###### Comet use py36, SC3 use yingluR4.1
    loadline = "module purge\nmodule load anaconda\nconda activate " + environment +"\n"
   # if environment =="":
   #     loadline = "module load anaconda\nconda activate gene_select\n"
   # elif environment == "sc3":
   #     loadline = "module load anaconda\nconda activate yingluR4.1\n"
   # else:
   #     loadline = "module load anaconda\nconda activate py36\n"
    timeline_start = "start=`date +%s`\n"
    comline = command
    timeline_end = "\nend=`date +%s`\nruntime=$((end-start))\necho $runtime second >> time_log"

    # Writing a string to file
    file.write(headline)
    file.write(loadline)
    file.write(timeline_start)
    file.write(comline)
    file.write(timeline_end)
    # Closing file
    file.close()

    sbatch_command = "sbatch --parsable sub.sh"
    sbatch_response = subprocess.getoutput(sbatch_command)
    job_id = sbatch_response.split(" ")[-1].strip()

    return job_id



# step1: preprocess
def preprocess(data_dir, work_dir, nvariable=2000, known_marker = False, keep_marker = False, Cluster=True, nPCA=10, max_RNA=2000, min_RNA=200, max_mt=5, partition="yaolab,batch", environment="gene_select"):

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
        + str(known_marker)
        + " "
        + str(keep_marker)
        + " "
        + str(max_RNA)
        + " "
        + str(min_RNA)
        + " "
        + str(max_mt)
        + " "
        + str(Cluster)
        + " "
        + str(nPCA)
        + " > stat_preprocess"
    )
    job_id = sbatch(job_name="preprocess", command=command, work_dir=work_dir, partition=partition, environment=environment)
    print(job_id)

    return job_id



# step 2: Selection


def selection(work_dir, data_dir="", method="de", n_marker=10, dep="", partition="yaolab,batch", environment="gene_select", **kwarg):

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

    job_id = options[method](work_dir=work_dir, data_dir=data_dir, n_marker=n_marker, dep=dep, partition=partition, environment=environment,  **kwarg)
    return job_id



# step 3: Evaluation

# supervised evaluation
def evaluation(work_dir, nPCA=10, truncat_n = 0, dep="", partition="yaolab,batch", environment="gene_select"):

    # cluster again
    command = "Rscript " + code_dir + "re-cluster.r " + work_dir  + " " + str(nPCA) + " " + str(truncat_n)
    # command = "python " + code_dir + "test.py"
    job_id_cluster = sbatch(job_name="re-cluster", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    print(job_id_cluster)

    # compare y_predict with y_true
    command = "python " + code_dir + "evaluation.py " + work_dir + " > stat_evaluation"

    job_id = sbatch(job_name="evaluation", command=command, work_dir=work_dir, dep=str(job_id_cluster), partition=partition, environment=environment)
    print(job_id)
    return



# seurat differential analyses method
def diff_express(work_dir, data_dir="", n_marker=10, dep="", partition="yaolab,batch", environment="gene_select"):
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
    job_id = sbatch(job_name="selection", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    return job_id



# SC3 method
# use envrionment yingluR4.1
def SC3_diff(work_dir,data_dir="", n_marker=10, dep="", partition="yaolab,batch", environment="gene_select"):
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
    job_id = sbatch(job_name="selection", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    return job_id



# scGenefit method
def scGenefit(
    work_dir, data_dir="",input_format="10X", n_marker=10, method="centers", epsilon= 1, redundancy=0.25, dep=""
, partition="yaolab,batch", environment="gene_select"):
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
    job_id = sbatch(job_name="selection", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    return  job_id



# Comet method
# use envrionment py36
def Comet(
    work_dir, data_dir="", n_marker=10, vis=False, if_pair=1, dep="", others=""
, partition="yaolab,batch", environment="gene_select"):
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
        tabmarker_file = data_dir + "/tabmarker.txt"
        tabcluster_file = data_dir + "/tabcluster.txt"
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
    print(command)
    job_id_comet = sbatch(job_name="selection",time=24, command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)
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

    job_id = sbatch(job_name="comet_result", time=8, command=command, work_dir=work_dir, dep=str(job_id_comet), partition=partition, environment=environment)
    
    return job_id



# SCmarker method
def SCmarker(work_dir, data_dir="", n_marker=10, k=100, n=10, dep="",partition="yaolab,batch", environment="gene_select"):

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
    job_id = sbatch(job_name="selection", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    return job_id


# COSG method
def COSGmarker(work_dir, data_dir="",n_marker=10, mu=1, dep="", partition="yaolab,batch", environment="gene_select"):
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
    job_id = sbatch(job_name="selection", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    return job_id

# FEAST method
def FEASTmarker(work_dir,data_dir="", n_marker=10, dep="", partition="yaolab,batch", environment="gene_select"):
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
    job_id = sbatch(job_name="selection", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    return job_id


# high vaiable genes method
def high_variable(work_dir, data_dir="", n_marker=10, dep="", partition="yaolab,batch", environment="gene_select"):

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
    job_id = sbatch(job_name="selection", command=command, work_dir=work_dir, dep=dep, partition=partition, environment=environment)

    return job_id

