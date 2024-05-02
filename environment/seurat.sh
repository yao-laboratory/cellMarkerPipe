#!/bin/bash

cp environment/seurat.yaml ./
conda env create -f seurat.yaml
conda activate seurat_env
pip install -e .
