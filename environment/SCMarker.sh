#!/bin/bash

cp environment/SCMarker.yaml ./
conda env create -f SCMarker.yaml
source activate SCMarker_env
# Run the R commands to install the SCMarker package
R -e "devtools::install_github('KChen-lab/SCMarker')"
pip install -e .
