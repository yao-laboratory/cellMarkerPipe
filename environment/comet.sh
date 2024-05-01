#!/bin/bash

# cp environment/comet.yaml ./
# conda env create -f comet.yaml
# conda activate COMET_env
# pip install -e .

cp comet.yaml ../
cd ..
conda env create -f comet.yaml
source activate COMET_env
pip install -e .
