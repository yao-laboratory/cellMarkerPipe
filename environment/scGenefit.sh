#!/bin/bash

# cp environment/scGenefit.yaml ./
# conda env create -f scGenefit.yaml
# conda activate scGenefit_env
# pip install -e .

cp scGenefit.yaml ../
cd ..
export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
conda env create -f scGenefit.yaml
source activate scGenefit_env
pip install scGeneFit
pip install -e .
