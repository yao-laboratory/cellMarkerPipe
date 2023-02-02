#!/bin/bash

conda create -n gene_select python=3.9
conda install -c conda-forge r=4.1 -n gene_select
conda install -c conda-forge r-seurat -n gene_select
conda install -c bioconda bioconductor-dropletutils -n gene_select
conda activate gene_select
pip install -e .
