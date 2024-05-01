#!/bin/bash

# cp environment/SC3.yaml ./
# conda env create -f SC3.yaml
# conda activate sc3_env
# pip install -e .

cp SC3.yaml ../
cd ..
conda env create -f SC3.yaml
source activate sc3_env

# Install R packages
R --vanilla << EOF
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SC3")
EOF
pip install -e .
