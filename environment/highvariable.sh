#!/bin/bash

cp environment/highvariable.yaml ./
conda env create -f highvariable.yaml
conda activate highvariable_env
pip install -e .
