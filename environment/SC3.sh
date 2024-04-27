#!/bin/bash

cp environment/SC3.yaml ./
conda env create -f SC3.yaml
conda activate sc3_env
pip install -e .
