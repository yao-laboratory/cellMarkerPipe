#!/bin/bash

cp environment/seurat.yaml ./
conda env create -f seurat.yaml
source activate seurat_env
pip install -e .
