# cellMarkerPipe
[![DOI:10.1038/s41598-024-63492-z](https://zenodo.org/badge/DOI/10.1038/s41598-024-63492-z.svg)](https://doi.org/10.1038/s41598-024-63492-z)
[![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.1038/s41598-024-63492-z)](https://doi.org/10.1038/s41598-024-63492-z)
[![GitHub stars](https://img.shields.io/github/stars/yao-laboratory/cellMarkerPipe.svg?style=social&label=Star&maxAge=2592000)](https://GitHub.com/yao-laboratory/cellMarkerPipe/stargazers/)

<img src="./cellMarkerPipe.png" width="200" height="250"/>

## 1. Overview
CellMarkerPip is a pipeline designed to include recent popular tools of identifying marker genes using scRNAseq data and evaluate the effects of selected markers. Now the tools that have been tested and included in CellMarkerPip are:

(1) FindAllMarkers by Seurat ([Paper](https://doi.org/10.1016/j.cell.2021.04.048))  
(2) scGenefit ([Paper](https://doi.org/10.1038/s41467-021-21453-4))  
(3) SCMarker ([Paper](https://doi.org/10.1371/journal.pcbi.1007445))  
(4) SC3 ([Paper](http://dx.doi.org/10.1038/nmeth.4236))  
(5) COMET ([Paper](https://doi.org/10.15252/msb.20199005))  
(6) COSG ([Paper](https://doi.org/10.1093/bib/bbab579))  

## 2. Installation
This pipeline supports different gene-selection methods writen by different language (python & R). Install the running environment for each of the above tools are critical.
##### Step 1: Download cellMarkerPipe
Please download the package directly from github (git clone). (or You may download our entire package as a zip file)
```shell
git clone https://github.com/yao-laboratory/cellMarkerPipe.git
```
##### Step 2: One-stop environment installation using bash (.sh) file
Go to the folder of cellMarkerPipe
Run bash file in environment/*.sh to install the method environment that you want to use.
For example, if you want to use the "FindAllMarkers method" or any marker gene selection algorithms provided by Seurat, you should use `seurat.sh` to set the envrionment. It will install all nessesary dependencies in the bash file.

```shell
bash environment/seurat.sh
```
##### (Alternative) Step 2: Step-by-step environment installation
This will be useful if running above bash file is difficult to accomplish in your system, or you need more freedom in specify environment name, or you need some flexibility.
Step-by-step installation can be easily followed up and done by looking into the steps inside environment/*.sh file

We use Seurat method as an example.

(1)
After you make sure anaconda is available, you should run command to install environment from yaml file and activate environment.
``` shell
cp environment/seurat.yaml ./
conda env create -f seurat.yaml
```
Then you should have a customized conda envrionment named seurat_env with cellMarkerPipe installed. To use cellMarkerPipe, you need to activate this envrionment.
``` shell
conda activate seurat_env
```
If you experience hard time in conda installation, mamba is a good subsitution.
``` shell
conda create -n my_env
conda activate my_env
conda install -c conda-forge mamba # if you haven't got mamba in your system
mamba env update -n my_env --file seurat.yaml
```
(2)
You may need to install additional packages for some method.

As for what additional packages are required excpet the ones listed in .yaml file, please check the .sh files corresponding to each method under the envrionment folder.

(3)
After you successfully install the method envrionment seurat_env, you can use pip to register this package 'cellMarkerPipe' to your own python path so that you can run commands in any location.
``` shell
pip install -e .
```
##### Step 3: Simple test for the installation
The enviroment and packages should be created and installed from above steps. And you can check whether the installation is successful by running cellMarkerPipe in command line from any directory. 
``` shell
cellMarkerPipe --version
```
If the package is successfully installed, the screen will show you the version of this package as
```
0.0.0
```

## 3. Tutorial
#### Input data
(required inputs)
`DATADIR` is a user defined data directory.

The neccessary input data is the cell matrix data in 10x format (`matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`) under `DATADIR`. These files can be prepared by 10x pipeline. 

If you have other data types such as Seurat data object, currently we provide a jupyternotebook example to convert from Seurat data object to 10x format. (please see testsuit/tutorial/ex2_start_from_seurat.ipynb)

(optional inputs)
You can choose to provide clustering(group) information of the cells if you have already clustered the cells using other tools or pipelines. This file needs to be provided in file named `groups.csv` under `DATADIR`. The `groups.csv` needs to contain 2 columns seperated by `","` the first column is the cell barcodes (same in `barcodes.tsv` ) and the second column is the cell cluster ID. 

If you have known marker genes for cell clusters, you may provide this file named as `Known_marker.csv` inside the data directory `DATADIR`

An example dataset is provided in the folder `data/Zeisel/10x/` inside this package.

Then you can run the pipeline as example below:

#### Key output data
The main result of this pipepline is:
(1) `marker_gene_per_group.csv` in marker folder includes marker genes for each group.
（2） `result.csv` in evaluation folder includes evaluation index of ari,jaccard,purity,nmi,fmi,f1_macro,f1_micro after the recluster is proceeded.
（3） `precision_recall.csv` in evaluation folder includes selected marker genes per group which are also shown in `Known_marker.csv` as listed in the first row with rowname "found_markers". The second row and the third row are precision and recall for each group. The last two rows include the total precision and recall.

#### Run cellMarkerPipe in command-line mode
##### Step 0: cellMarkerPipe Commands Overview
This pipepline has 3 main steps: preprocess, selection and evaluation. 

``` bash
cellMarkerPipe -h
```
If the pacakge is successfully installed, you should find this output on your screen
```
usage: cellMarkerPipe [-h] [--version] {preprocess,selection,evaluation} ...

Find marker genes for single cell datasets

positional arguments:
  {preprocess,selection,evaluation}
                        help for subcommand: preprocess, selection, evaluation
    preprocess          Preprocess the 10x data
    selection           Select marker genes
    evaluation          Evaluate selected marker genes

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```
##### Step 1: Preperation
This step is to use single cell RNA seq data to do preparations for marker gene selection, such as pre-processing and filtering by RNA quality, mitochondra genes, dimension reduction, and cell clustering etc.
For the first step `preperation`, you can check how to set up the parameters by run:

``` shell
cellMarkerPipe preprocess -h
```
The parameters that you can set up are shown in output:

```
usage: cellMarkerPipe preprocess [-h] -wd WORKDIR -10xd DATADIR [-nvb NVARIABLE] [-maR MAXRNA] [-miR MINRNA] [-mam MAXMT]
                                 [-np NPCA] [-res RESOLUTION] [-alg ALGORITHM] [--cluster] [--no-cluster] [--know-marker]
                                 [--no-know-marker] [--keep-known-marker] [--no-keep-known-marker]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -10xd DATADIR, --10xdir DATADIR
                        10x data directory
  -nvb NVARIABLE, --nvariable NVARIABLE
                        The number of highly variable genes selected for later selection
  -maR MAXRNA, --maxRNA MAXRNA
                        max number of RNA for each cell
  -miR MINRNA, --minRNA MINRNA
                        min number of RNA for each cell
  -mam MAXMT, --maxmt MAXMT
                        max number of MT genes for each cell
  -np NPCA, --nPCA NPCA
                        The number of PCA chosen for re-cluster
  -res RESOLUTION, --resolution RESOLUTION
                        The solution value used in FindClusters for re-cluster
  -alg ALGORITHM, --algorithm ALGORITHM
                        The algorithm chosen in FindClusters for re-cluster
  --cluster             Do clustering
  --no-cluster          Do not do clustering, you already have cluster information
  --know-marker         Have a file named Known_marker.csv
  --no-know-marker      Do not have a file named Known_marker.csv
  --keep-known-marker   Keep the known markers during screening
  --no-keep-known-marker
                        Do not keep the known markers during screening
```
Among these parameters, `-wd WORKDIR` (where you want to save outputs) and `-10xd DATADIR` (where you have inputs) are two parameters that are necessary to be provided. 

In this tutorial, we use the `testsuit/test` folder under the package `cellMarkerPipe` as the `WORKDIR` to have a testrun. There is an simple Zeisel 10x dataset provided under the folder `data/Zeisel/10x` in the package as the `DATADIR`


``` bash
cd testsuit/test
#The relative path is used in this example. We also strongly suggest to provide absolute paths.
cellMarkerPipe preprocess -wd ./ -10xd ../../data/Zeisel/10x
#
```


After this step is finished, a folder named `Data` and a file named `stat_preprocess` will be generated under the `WORKDIR`. The file `stat_preprocess` includes the standard ouput of the program, while `Data` folder include the subset of count matrix of high variable genes in 10X and csv format, which are required for next `Selection` step. 

There are some optional parameters provided to users to customize the data screening and cluster process, where we followed the tutorial of Seurat (https://satijalab.org/seurat/articles/pbmc3k_tutorial) to conduct:
1. coarsely filter cells that have unique feature counts over `MAXRNA` or less than `MINRNA`. We filter cells that have >`MAXMT` mitochondrial counts. A meta-ouput of distribution_of_features_counts.png provide you with a reference to choose appropriate values for your dataset.
2. data normalization,
3. identify high-variable genes, and only use top `NVARIABLE` genes for following calculation
4. scale data
5. reduce dimenion with PCA techinique. `NPCA` principle components are used for following cluster
6. apply a graph-based clustering approach. You can customize the parameters `RESOLUTION` and `ALGORITHM` of Seurat method FindClusters used.

In above process, step 5&6 are not neccessary to do if you already provide cell clusters in the file `groups.csv` manually or by other clustering tools. Then make sure `--no-cluster` is set. Otherwise you can use `--cluster` to allow cellMarkerPipe run clustering for you. By default, cellMarkerPipe will do clustering automatically.

CellMarkerPipe can be used to compare the selected marker genes with well-known marker genes. So we also provide the users a choice whether keep all well-known marker genes in case they are not in high variable gene list. To keep all well-known genes, you need to provide a file named as 'Known_marker.csv' under the `DATADIR` and list the well-known genes for each cluster in each row. The 1st column is the cluster ID or label, and the second column is the known marker genes separated by comma. By default, the pipeline will not keep the well-known genes, you need to mannually adding parameters `--know-marker` and `--keep-known-marker` to keep them in high variable gene list and for evaluation step.

##### Step 2: Select Marker Genes
This step is to select marker genes by different methodology. Different methods have been wrapped by python or R codes in cellMarkerPipe.
To run `selection` step, you can check how to set up the parameter by run,

``` shell
cellMarkerPipe selection -h
```
Then you should see all parameters that you can set up from the output:

```
usage: cellMarkerPipe selection [-h] -wd WORKDIR -m METHOD [-n NMARKER]
optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -m METHOD, --method METHOD
                        Method used for selection
  -n NMARKER, --n-marker NMARKER
                        Number of marker genes to select
```
Here, `WORKDIR` should keep the same as `preperation` step. Then you need to choose the `METHOD` you want to use. Here are the mapping between the method and its shortname:
| Method | FindAllMarkers    | scGenefit    | COSG | SC3 |COMET|SCMarker|High Variable|FEAST|
| :---:   | :---: | :---: |:---: |:---: |:---: |:---: |:---: |:---: |
| Shortname | de   | scG   | cos | sc3 | Com | SC | hv | fst |

Another parameter `NMARKER` decides the number of marker genes assigned to each cluster. The default value is 10. 

For example, we using method `de` (FindAllMarkers in Seurat) in this example command.
``` bash
cellMarkerPipe selection -wd ./ -m de
```
Using the example data, this step takes about a few minutes depending on which method you use. After this step is finished, you may find a `marker` folder under your `WORKDIR`, which include the results of selection. `marker_gene_per_group.csv` includes the selected marker genes for each group, which is required if you want to do the `Evaluation` step. The other files are meta-data for your reference. For example `stat_selection` is provided for you under `WORKDIR` to debug.

##### Step 3: Evaluation
This step is to evaluation the quality of the marker genes selected by above methods and provide metrics. 
To understand the parameters, you can run
```
cellMarkerPipe evaluation -h
```
Then you find all explanations of the parameters by output:
```
usage: cellMarkerPipe evaluation [-h] -wd WORKDIR [-np NPCA] [-res RESOLUTION]
                                 [-alg ALGORITHM] [--know-marker]
                                 [-10xd DATADIR]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -np NPCA, --nPCA NPCA
                        The number of PCA chosen for re-cluster
  -res RESOLUTION, --resolution RESOLUTION
                        The solution value used in FindClusters for re-cluster
  -alg ALGORITHM, --algorithm ALGORITHM
                        The algorithm chosen in FindClusters for re-cluster
  --know-marker         Have a file named Known_marker.csv
  -10xd DATADIR, --10xdir DATADIR
                        10x data directory
```
Here you need to provide  `WORKDIR`, which should be the same as above steps. In the process of `evaluation`, the program will redo the cluster with only the selected marker genes from the `selection` step to check the quality of marker genes, which is called re-clustering. Same as `preperation` step, you can customize the clustering by setting up the nearest neighbor graph (see Seurat paper or document) by `nPCA`, `RESOLUTION` and `ALGORITHM` of FindClusters method in seurat. 

You can run example command below.
```bash
cellMarkerPipe evaluation -wd ./ -np 10
```

The outputs can be found under folder `evaluation` with a filename `result.csv`. This step also create a `re-cluster` folder under the `WORKDIR` which is the output of the re-clustering process.
The final evaluation report includes re-clustering scores such as the Adjusted Rand Index (ARI),Jaccard index, purity, normalized mutual information (NMI), Fowlkes-Mallows Index (FMI) and Marco-F1 and Micro-F1 scores. If the user provides known marker genes, additional score report will be provided for precision and recall values for each cell type and overall dataset comparing to known marker genes (`precision_recall.csv`). To switch on the calcuation of precision and recall, you need to use parameter --know-marker to switch on the calculation and provide the progam the directory where you store the known marker genes ``Known_marker.csv`` after the parameter -10xd. The `Known_marker.csv` is the same file as we discussed in section `Preparation`.

#### (alternative) Run cellMarkerPipe as a library in your own python codes or in Jupyternotebook
For developers, we also provide you a method to use cellMarkerPipe as a python library. 

##### Step 0: Import pipeline and define data_dir and work_dir
``` python
from pipeline import block
import os
# input data directory: 10X. Same parameter as "DATADIR" in command line mode
data_dir = "YOUR_DATA_DIR"
# output data directory. Same parameter as "WORKDIR" in command line mode
work_dir = "YOUR_WORK_DIR"
if not os.path.isdir(work_dir):
    os.mkdir(work_dir)
os.chdir(work_dir)
```
##### Step 1: Preperation
``` python
block.preprocess(data_dir=data_dir, work_dir=work_dir, nvariable=2000, Cluster=False, max_RNA = 2500, max_mt = 5)
```
Here `navriable`, `max_RNA` and `max_mt` are the same parameters corresponding to `NVARIABLE`, `MAXRNA` and `MAXMT` inseperately  command line. "Cluster" decides whether you need to do cluster in this step (same as `--cluster` and `--no-cluster`). If you do not provide `group.csv`, then this parameter has to be turn to "TRUE", which will execute Seurat clustering cells. 

##### Step 2: Select Marker Genes
``` python
block.selection(work_dir=work_dir, method="de")
```
Here `method` is the same parameter as `METHOD` in command line 
##### Step 3: Evaluation
``` python
block.evaluation(work_dir, nPCA=10)
```
`nPCA` is the same parameter as `NPCA` in command line
The above steps have been included into the test file under folder 'notebook'. The test data are under folder 'data'.

## 4. Additional Topics

#### 4.1 Support for Seurat Object as input.

Officially, we only support 10x sparse matrix format in the pipeline. In order to use your seurat object in this pipeline, we suggest you to convert it to the 10x format.
We have provided an example data (R object) in the folder `data/Zeisel/seurat/seurat.rds`
Please follow the jupyternotebook to do the conversion. `testsuit/tutorial/ex2_start_from_seurat.ipynb`

#### 4.2 Inject known marker genes as ground truth for evaluation.

If you have known marker genes for cell clusters, you may provide this file named as `Known_marker.csv` inside the data directory `DATADIR`. This is already explained in step 1. By default, the pipeline will not consider the well-known genes; you need to mannually adding parameters `--know-marker` and `--keep-known-marker` to consider them. 

#### 4.3 Inject marker genes selected from other computational tools to evaluate.

First you can run the entire pipeline with your cell matrix data by any of our provided method, without allowing our pipeline to do clustering for you in the preparation. In the second step, you will see a selected marker gene file `marker_gene_per_group.csv`. Simply modify this file to include marker genes from your computational tool, and then run the evaluation step.

#### 4.4 Integrate a new tool completely into the pipeline (as a developer).
Developers are welcome to integrate the gene selection codes writen in Python or R into this pipeline. To achieve that, 

Here are the procedures to help you achieve that.
##### Step 1: prepare a Python/R script to include the gene selection codes
###### Example standard input and output for Python: YOUR_METHOD.py
``` Python
#!/usr/bin/env python
# coding: utf-8

# import libarary
import ANY_LIBRARY_YOU_NEED
########################### input #############################################################################
# The standard selection method requires the script must accept 3 arguments:
# work_dir: the work directory
# n_marker: the number of marker genes selected for each group

# input variables
work_dir = sys.argv[1]
print("Work directory is:")
print(work_dir)

n_marker = sys.argv[2]
n_marker = int(n_marker)
print("The number of markers is:")
print(n_marker)

# read gene*cell matrix into a sc.adata object
# it is recommended to import 10x data because the size is small, so the import process takes little time
# Actually, there are normalized and scaled matrix provided in work.dir/data/ as the meta-data file
# You can also import those files if you have a small dataset.

data_path=work_dir + "/data/10x"
    adata = sc.read_10x_mtx(
    data_path,  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=False)
genes = adata.var_names

# import cluster lable file
cluster_dir = work_dir + "/data/cluster_labels.csv"
clusters = pd.read_csv(tabcluster_file, sep = '\t', header=None)
clusters_np = clusters.iloc[:,1].values
Y_t, indices = np.unique(clusters_np, return_inverse=True)

#############################################################################################################
############################### body ########################################################################
# Here you customize the file by add the code neccesary to deal with the matrix and group labels to select
# marker genes. Finally, you should provide a dataframe marker_per_group to save the marker genes for each group
marker_per_group = pd.DataFrame(index=Y_t, columns=["Marker"])
marker_per_group.index.name = "Cluster"
############################################################################################################
############################### output #####################################################################
select_dir = work_dir + "/marker"
result_per_group_dir = select_dir + "/marker_gene_per_group.csv"
marker_per_group.to_csv(result_per_group_dir, header = True, sep=',')
############################################################################################################
```

###### Example standard input and output for R: YOUR_METHOD.R
``` R
# import libarary
library(ANY_LIBRARY_YOU_NEED)

########################### input #############################################################################
# The standard selection method requires the script must accept 3 arguments:
# work.dir: the work directory
# n.marker: the number of marker genes selected for each group

args <- commandArgs(trailingOnly = TRUE)

## start a project dir
work.dir = args[1]
n.marker = args[2]

print("The work directory is:")
print(work.dir)

print("The number of marker chosen for each group:")
n.marker <- as.integer(n.marker)
print(n.marker)

# read gene*cell matrix
# it is recommended to import 10x data because the size is small, so the import process takes little time
# Actually, there are normalized and scaled matrix provided in work.dir/data/ as the meta-data file
# You can also import those files if you have a small dataset.
data.dir <- file.path(work.dir, "/data/10x")
pbmc.data <- Read10X(data.dir = data.dir)

# read group lables for each cell
cluster.dir <- file.path(work.dir, "data/cluster_labels.csv")
cluster <- read.table(
  file = cluster.dir,
    sep="\t",
  as.is = TRUE
)

###############################################################################################################

########################### body ##############################################################################
# Here you customize the file by add the code neccesary to deal with the matrix and group labels to select
# marker genes. Finally, you should provide a dataframe marker.per.group to save the marker genes for each group
marker.per.group <- data.frame(Cluster=character(0), Marker=character(0))

###############################################################################################################

########################### output ############################################################################
# Save results
marker.dir <- file.path(work.dir, "marker")
if (!file.exists(marker.dir)) {
dir.create(marker.dir)
}
write.csv(marker.per.group, file.path(marker.dir, "marker_gene_per_group.csv"),row.names = FALSE)
###############################################################################################################
```

##### Step 2: Add the name of the srcipt into the list
In the file `cellMarkerPipe/pipeline/block.py` you can find the `selection` method. You need to add your own script's name into the dictionary using abbreviation as the key and the name fo the python/R script as the value.
```
def selection(work_dir, method="de", n_marker=10,  **kwarg):

    # map the inputs to the function blocks
    options = {
        "de": diff_express,  # by group
        "scG": scGenefit,  # not by group
        "Com": Comet,  # by group, by neighbouring number
        "SC": SCmarker,  # not by group
        "cos": COSGmarker,  # by group
        "fst": FEASTmarker,  # not by group
        "sc3": SC3_diff, # by group
        "hv": high_variable, # by group
        "abbreviation"： YOUR_METHOD # your own python/R script name
    }

    options[method](work_dir=work_dir,  n_marker=n_marker,  **kwarg)
    return
```
##### Step 3: write a function named YOUR_METHOD in block.py
Then you need to write a function named method to transfer parameters required by `YOUR_METHOD.R/YOUR_METHOD.py` in `block.py`. Here is an example:
```python
# YOUR_METHOD
def YOUR_METHOD(work_dir, n_marker=10):
    # n is the number of marker for each cluster
    command = (
        "Rscript "
        + code_dir
        + "YOUR_METHOD.r "
        + work_dir
        + " "
        + str(n_marker)
        + " > stat_selection"
    )
    os.system(command)
    return 
```


## Citation

If you use this work, please cite:

```bibtex
@article{jia2024cellmarkerpipe,
  title={CellMarkerPipe: cell marker identification and evaluation pipeline in single cell transcriptomes},
  author={Jia, Y., Ma, P. and Yao, Q.},
  journal={Scientific Reports},
  volume={14},
  pages={13151},
  year={2024},
  publisher={Nature Publishing Group},
  doi={10.1038/s41598-024-63492-z},
  url={https://doi.org/10.1038/s41598-024-63492-z}
}
