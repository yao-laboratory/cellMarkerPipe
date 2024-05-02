#!/bin/bash

cp environment/highvariable.yaml ./
conda env create -f highvariable.yaml
source activate highvariable_env
# Run the R commands to install the SCMarker package
R -e "devtools::install_github('KChen-lab/SCMarker')"
pip install -e .
