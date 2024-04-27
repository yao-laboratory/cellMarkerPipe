#!/bin/bash

cp environment/cosg.yaml ./
conda env create -f cosg.yaml
conda activate COSGR_env
pip install -e .
