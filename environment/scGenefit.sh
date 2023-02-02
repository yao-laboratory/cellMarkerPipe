#!/bin/bash

conda create -n select2 python=3.9
conda install -c conda-forge r=4.1 -n select2
conda install -c conda-forge r-seurat -n select2
conda install -c bioconda bioconductor-dropletutils -n select2
conda activate select2
pip install scanpy 

# may need to update pip
# pip install update

pip install scGenefit
pip install -e . 
