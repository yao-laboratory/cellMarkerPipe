#!/bin/bash

# cp environment/cosg.yaml ./
# conda env create -f cosg.yaml
# conda activate COSGR_env
# pip install -e .

cp cosg.yaml ../
cd ..
conda env create -f cosg.yaml
source activate COSGR_env

# Install R packages
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')"
Rscript -e "remotes::install_github('genecell/COSGR')"
pip install -e .
