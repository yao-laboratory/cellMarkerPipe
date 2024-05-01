#!/bin/bash

# cp environment/seurat.yaml ./
# conda env create -f seurat.yaml
# conda activate seurat_env
# pip install -e .

cp seurat.yaml ../
cd ..
conda env create -f seurat.yaml
source activate seurat_env
pip install -e .
