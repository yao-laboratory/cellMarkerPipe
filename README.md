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
The enviroment and software should be created and installed after this step. And you can check whether the installation is successful by running cellMarkerPipe in command line. 
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
python cellMarkerPipe.py -h
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
To run `preperation` step, you can use command.
``` bash
python cellMarkerPipe.py preprocess -wd ./ -10xd data/Zeisel/10x
```
```
usage: cellMarkerPipe preprocess [-h] [-wd WORKDIR] [-10xd DATADIR] [-nvb NVARIABLE] [--cluster] [--no-cluster] [-maR MAXRNA]
                                 [-mam MAXMT]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -10xd DATADIR, --10xdir DATADIR
                        10x data directory
  -nvb NVARIABLE, --nvariable NVARIABLE
                        The number of highly variable genes selected for later selection
  --cluster             Do cluster
  --no-cluster          Do not do cluster
  -maR MAXRNA, --maxRNA MAXRNA
                        max number of RNA for each cell
  -mam MAXMT, --maxmt MAXMT
                        max number of MT genes for each cell
```

In the example above, the relative path is used for `WORKDIR` and `DATADIR`.  However, it is strongly suggested to provide parameters `WORKDIR`, `DATADIR` with an absolute path. After this step is finished, a folder named `Data` and a file named `stat_preprocess` will be generated under the `WORKDIR`. The file `stat_preprocess` includes the standard ouput of the program. While `Data` folder include the subset of count matrix of high variable genes in 10X and csv format, which are required for next `Selection` step. The number of high variable genes you want to study with in the next step can be set up using `NVARIABLE` parameter. `--cluster` or `--no-cluster` decides whether you need the programe to do cluster process firstly, which depends on whether you provide input file of `groups.csv`.  A meta-ouput of distribution_of_features_counts.png provide you with a reference to coarsely filter genes coarsely by choosing appropriate truncate of max number of RNA for each cell (`MAXRNA`) and max number of MT genes for each cell(`MAXMT`), which are also important to keep the informative genes while clean genes and cell at the tail of the distribution to eliminating noise to the dataset. Using the example data, this step takes less than 10 seconds. 

##### Step 2: Select Marker Genes
To run `selection` step, you can use command, we using method `de` in this example command.
``` bash
python cellMarkerPipe.py selection -wd ./ -10xd /.../cellMarkerPipe/data/Zeisel/10x -m de
```
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

Using the example data, this step takes about a few minutes depending on which method you use. After this step is finished, you may find a `marker` folder under your `WORKDIR`, which include the results of selection. `marker_gene_per_group.csv` includes the selected marker genes for each group, which is required if you want to do the `Evaluation` step. The other files are meta-data for you reference.
##### Step 3: Evaluation
We also provide users a unsurpervised method to evalute the selected markers by calculating indexs which evalute how these marker genes can seperate the cell, including ARI et. al. You can use example command below.
```bash
python cellMarkerPipe.py evaluation -wd ./ -np 10
```
```
usage: cellMarkerPipe evaluation [-h] [-wd WORKDIR] [-np NPCA]

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKDIR, --workdir WORKDIR
                        Working directory
  -np NPCA, --nPCA NPCA
                        The number of PCA chosen for re-cluster
```
This evaluation needs to redo the cluster process using only the selected markers. So you can set up `nPCA` to optimize this cluster process. `WORKDIR` should keep the same.  Using the example data, this step takes about less than 10 seconds. After you finished this step, the calculated index can be find under folder `evaluation` with a filename `result.csv`. the program also create a `re-cluster` folder under the `WORKDIR` which is the output of the cluster process.
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
