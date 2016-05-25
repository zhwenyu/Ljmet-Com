#!/bin/bash

echo "SUBMITTING LJMET -- PRD Data"

cp TTtrilepton_Data_cfg.py ljmet_cfg.py

python -u condor_submitargs_data.py nominal

echo "DONE"