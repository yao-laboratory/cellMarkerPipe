# cellMarkerPipe

<img src="./cellMarkerPipe.png" width="200" height="250"/>

### Overview
CellMarkerPip is a pipeline designed to include recent popular tools of identifying marker genes using scRNAseq data and evaluate the effects of selected markers. Now the tools that have been built in CellMarkerPip are:

(1) FindAllMarkers by Seurat ([Paper](https://doi.org/10.1016/j.cell.2021.04.048))  
(2) scGenefit ([Paper](https://doi.org/10.1038/s41467-021-21453-4))  
(3) SCMarker ([Paper](https://doi.org/10.1371/journal.pcbi.1007445))  
(4) SC3 ([Paper](http://dx.doi.org/10.1038/nmeth.4236))  
(5) COMET ([Paper](https://doi.org/10.15252/msb.20199005))  
(6) COSG ([Paper](https://doi.org/10.1093/bib/bbab579))  

### Installation
This pipeline supports different gene-selection method writen by different language (python & R). After download the package to your cluster, firstly use anaconda to set up envrionment, then install the package with pip install. These steps have been included into installing files under folder `envrionment`. Some methods have dependence conflict. You have to use the right one to set up environment. For example, if you want to use seurat method, then you should run command and activate environment.
``` shell
git clone https://github.com/yao-laboratory/cellMarkerPipe.git
cd cellMarkerPipe
cp environment/seurat.yaml ./
conda env create -f seurat.yaml
```
Then you should have a customized conda envrionment named seurat_env with cellMarkerPipe installed. To use cellMarkerPipe, you firstly need to activate this envrionment.
``` shell
conda activate seurat_env
```
If user experiences hard time in conda installation, mamba is a good subsitution.
``` shell
conda create -n my_env
conda activate my_env
conda install -c conda-forge mamba # if you haven't got mamba in your system
mamba env update -n my_env --file seurat.yaml
```
Under the envrionment seurat_env, you can use pip to install this package under the folder cellMarkerPipe.
``` shell
pip install -e .
```
The enviroment and software should be created and installed after this step. And you can check whether the installation is successful by running cellMarkerPipe in command line at any directory. 
``` shell
cellMarkerPipe --version
```
If the package is successfully installed, the screen will show you the version of this package as
```
0.0.0
```

### Tutorial
#### Input
The neccessary input file is the counts matrix data in 10x format (`matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`) under `DATADIR`. You can choose to provide group information of the cells or not. If you want to use your own cell group, then it needs to be provided in file named `groups.csv` under `DATADIR`. The `groups.csv` needs to contain 2 columns seperated by `","` the first one is the cell barcodes same with `barcodes.tsv` and the second one is the cell name. An example dataset is provided in the folder `data/Zeisel/10x/` with the package.

Then you can run the pipeline as example below:

#### Using cellMarkerPipe in command-line mode
##### Step 0: cellMarkerPipe Overview
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
  --cluster             Do cluster
  --no-cluster          Do not do cluster
  --know-marker         Have a file named Known_marker.csv
  --no-know-marker      Do not have a file named Known_marker.csv
  --keep-known-marker   Keep the knwon markers during screening
  --no-keep-known-marker
                        Do not keep the known markers during screening
```
Among these parameters, `WORKDIR` and `DATADIR` are two parameters that are necessary to be provided with the path of your directories to save output and input 10x data, seperately. Here I use the `testsuit/test` folder under the package folder `cellMarkerPipe` as the `WORKDIR` to have a testrun. There is an simple Zeisel 10x dataset provided you under the folder `data/Zeisel/10x` kept in the package folder `cellMarkerPipe` for this testrun.

``` bash
cd testsuit/test
cellMarkerPipe preprocess -wd ./ -10xd ../../data/Zeisel/10x
```
In the example above, the relative path is used for `WORKDIR` and `DATADIR`.  However, it is strongly suggested to provide parameters `WORKDIR`, `DATADIR` with an absolute path. 

After this step is finished, a folder named `Data` and a file named `stat_preprocess` will be generated under the `WORKDIR`. The file `stat_preprocess` includes the standard ouput of the program. While `Data` folder include the subset of count matrix of high variable genes in 10X and csv format, which are required for next `Selection` step. 

There are some optional parameters provided to users to customize the data screening and cluster process, where we followed the tutorial of Seurat (https://satijalab.org/seurat/articles/pbmc3k_tutorial) to conduct:
1. coarsely filter cells that have unique feature counts over `MAXRNA` or less than `MINRNA`. We filter cells that have >`MAXMT` mitochondrial counts. A meta-ouput of distribution_of_features_counts.png provide you with a reference to choose appropriate values for your dataset.
2. data normalization,
3. identify high-variable genes, and only use top `NVARIABLE` genes for following calculation
4. scale data
5. reduce dimenion with PCA techinique. `NPCA` principle components are used for following cluster
6. apply a graph-based clustering approach. You can customize the parameters `RESOLUTION` and `ALGORITHM` of Seurat method FindClusters used.

In above process, step 5&6 are not neccessary to do if you already annote cluster labels to each cell. Then you need to prepare the file `groups.csv` which contains 2 columns seperated by `","` the first one is the cell barcodes same with `barcodes.tsv` and the second one is the cell name. We use `--cluster` or `--no-cluster` to switch whether cellMarkerPipe finds clusters for you. By default, the program does not find clusters automatically.

In our paper, we want to compare the selected marker genes with well-known marker genes. So we also provide the users a choice whether keep all well-known marker genes in case they are not top high variable genes and filtered by the pipeline. To keep all well-known genes, you need to provide a file named as 'Known_marker.csv' under the `DATADIR` and list the well-known genes for each cluster in each row. The 1st column is the cluster label, which can be anything and will be neglected later. By default, the pipeline will not keep the well-known genes, you need to mannually adding parameter --know-marker and --keep-known-marker to keep them for the next Selection step.

##### Step 2: Select Marker Genes
To run `selection` step, you can check how to set up the parameter by run,

``` shell
cellMarkerPipe selection -h
```
Then you should see all parameters that you can set up from the output:

```
usage: cellMarkerPipe selection [-h] [-wd WORKDIR] [-10xd DATADIR] [-m METHOD]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -10xd DATADIR, --10xdir DATADIR
                        10x data directory
  -m METHOD, --method METHOD
                        Method used for selection

```
Here, `WORKDIR` and `DATADIR` should keep the same as `preperation` step. Then you need to choose the method you want to use. Here are the mapping between the method and its shortname:
| Method | FindAllMarkers    | scGenefit    | COSG | SC3 |COMET|SCMarker|High Variable|FEAST|
| :---:   | :---: | :---: |:---: |:---: |:---: |:---: |:---: |:---: |
| Shortname | de   | scG   | cos | sc3 | Com | SC | hv | fst | 

For example, we using method `de` (FindAllMarkers in Seurat) in this example command.
``` bash
cellMarkerPipe selection -wd ./ -10xd ../../data/Zeisel/10x -m de
```
Using the example data, this step takes about a few minutes depending on which method you use. After this step is finished, you may find a `marker` folder under your `WORKDIR`, which include the results of selection. `marker_gene_per_group.csv` includes the selected marker genes for each group, which is required if you want to do the `Evaluation` step. The other files are meta-data for you reference. A standard output file `stat_selection` is provided for you under `WORKDIR` to debug.

##### Step 3: Evaluation
We also provide users a unsurpervised method to evalute the selected markers by calculating indexs which evalute how these marker genes can seperate the cell, including ARI et. al. To understand the parameters, you can run
```
cellMarkerPipe evaluation -h
```
Then you find all explanations of the parameters by output:
```
usage: cellMarkerPipe evaluation [-h] -wd WORKDIR [-np NPCA] [-res RESOLUTION] [-alg ALGORITHM]

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
```
Here you need to provide the code with  `WORKDIR`, which should keep the same as above steps. In the process of `evaluation`, the program will redo the cluster with only the selected marker genes from the `selection` step. As `preperation` step, you can customize the clustering by setting number of principle component used `nPCA`, `RESOLUTION` and `ALGORITHM` of FindClusters method in seurat. Then use the clustering result, the program can calculate the indexs like ARI et. al.

You can use example command below.
```bash
cellMarkerPipe evaluation -wd ./ -np 10
```

This evaluation needs to redo the cluster process using only the selected markers. So you can set up `nPCA` to optimize this cluster process.  Using the example data, this step takes about less than 10 seconds. After you finished this step, the calculated index can be find under folder `evaluation` with a filename `result.csv`. the program also create a `re-cluster` folder under the `WORKDIR` which is the output of the cluster process.
#### Using cellMarkerPipe as a library
For developers, we also provide you a method to use cellMarkerPipe as a library. The way how to use it is similar to the command line mode.

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
block.selection(work_dir=work_dir, data_dir ="", method="de")
```
Here `method` is the same parameter as `METHOD` in command line 
##### Step 3: Evaluation
``` python
block.evaluation(work_dir, nPCA=10)
```
`nPCA` is the same parameter as `NPCA` in command line
The above steps have been included into the test file under folder 'notebook'. The test data are under folder 'data'.

#### Integrate a new tool into the pipeline.
Developers are welcome to integrate the gene selection codes writen in Python or R into this pipeline. To achieve that, 

Here are the procedures to help you achieve that.
##### Step 1: prepare a Python/R script to include the gene selection codes
###### Example header and tail for Python

###### Example header and tail for R

##### Step 2: Add the name of the srcipt into the list
```
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
```
