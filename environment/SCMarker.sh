#!/bin/bash

cp environment/SCMarker.yaml ./
conda env create -f SCMarker.yaml
conda activate SCMarker_env
pip install -e .
