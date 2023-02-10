# cellMarkerPipe

### Overview
CellMarkerPip is a pipeline designed to include recent popular tools of identifying marker genes using scRNAseq data and evaluate the effects of selected markers. Now the tools that have been built in CellMarkerPip are:  
(1) FindAllMarkers by Seurat https://doi.org/10.1016/j.cell.2021.04.048  
(2) scGenefit https://doi.org/10.1038/s41467-021-21453-4  
(3) SCMarker https://doi.org/10.1371/journal.pcbi.1007445  
(4) SC3 http://dx.doi.org/10.1038/nmeth.4236  
(5) COMET https://doi.org/10.15252/msb.20199005  
(6) COSG https://doi.org/10.1093/bib/bbab579  

### Installation
This pipeline supports different gene-selection method writen by different language (python & R).  
After download the package to your cluster, firstly use anaconda to set up envrionment, then install the package with pip install.  
These steps have been included into installing files under folder "envrionment". Some methods have dependence conflict. You have to use the right one to set up environment.
For example, if you want to use seurat method, then you should copy the seurat.sh under the cellMarkerPipe directory and run
``` shell
bash seurat.sh
```

### Tutorial

##### Input
The neccessary input file is the counts matrix data in 10x format (`matrix.mtx`, `gene.tsv` and `barcode.tsv`) under 'data_dir'. You can choose to provide group information of the cells or not. If you want to use your own cell group, then it needs to be provided in file named `group.csv` under 'data_dir'. The "group.csv" needs to contain 2 columns seperated by ",": the first one is the cell barcodes same with `barcode.tsv` and the second one is the group name. Then you can run the pipeline as example below:

##### Import pipeline and define data_dir and work_dir
``` python
from pipeline import block
import os
# input data directory: 10X
data_dir = "YOUR_DATA_DIR"
# output data directory
work_dir = "YOUR_WORK_DIR"
if not os.path.isdir(work_dir):
    os.mkdir(work_dir)
os.chdir(work_dir)
```
##### Preperation
``` python
block.preprocess(data_dir=data_dir, work_dir=work_dir, nvariable=2000, Cluster=False, max_RNA = 2500, max_mt = 5)
```
The step is going to output a 'data' folder under your 'work_dir', which include meta-data needed for the next steps.  
Here "max_RNA" and "max_mt" are parameters to filter genes coarsely. You can choose propariate value by visualize the distribution of genes and cells, which can be found under output "data" folder. "nvariable" is the number of selected highly variable genes for selection in the next step. This parameter cut off the number of genes for selection method to work on. "Cluster" decides whether you need to do cluster in this step. If you do not provide `group.csv`, then this parameter has to be turn to "TRUE", which will execute Seurat clustering cells. Using the example data, this step takes less than 10 seconds. 

##### Select Marker Genes
``` python
block.selection(work_dir=work_dir, data_dir ="", method="de")
```
This step is going to ouput a 'marker' folder, in which you can find the marker genes for each group in file `marker_gene_per_group.csv`.
The parameter method decide which method you want to use to perform maker-gene selection. Using the example data, this step takes about a few minutes depending on which method you use.
Here are the mapping between the method and its shortname:
| Method | FindAllMarkers    | scGenefit    | SCMarker | SC3 |COMET|COSG|
| :---:   | :---: | :---: |:---: |:---: |:---: |:---: |
| Shortname | de   | scG   | SC | sc3 | Com | cos | 
 
##### Evaluation
If you want to evaluate the marker genes selected in the last step, we provided the unsupervised method to calculate indexs to evalute how these marker genes can seperate the cell, including ARI et. al. 
``` python
block.evaluation(work_dir, nPCA=10)
```
You can find out the calculated index under folder 'evaluation'.  Using the example data, this step takes about less than 10 seconds. 
The above steps have been included into the test file under folder 'notebook'. The test data are under folder 'data' while the corresponding output are under folder 'ouput'.
