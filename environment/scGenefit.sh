#!/bin/bash

cp environment/scGenefit.yaml ./
conda env create -f scGenefit.yaml
conda activate scGenefit_env
pip install -e .
