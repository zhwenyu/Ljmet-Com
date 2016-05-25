#!/bin/bash

# echo "SUBMITTING LJMET -- PRD Data"
# 
# cp TTtrilepton_Data_cfg.py ljmet_cfg.py
# 
# python -u condor_submitargs_data.py nominal

echo "SUBMITTING LJMET -- nominal"

cp TTtrilepton_MC_cfg.py ljmet_cfg.py

python -u condor_submitargs.py nominal

# echo "SUBMITTING LJMET -- JECup"
# 
# python -u condor_submitargs.py JECup
# 
# echo "SUBMITTING LJMET -- JECdown"
# 
# python -u condor_submitargs.py JECdown

echo "SUBMITTING LJMET -- JERup"

python -u condor_submitargs.py JERup

echo "SUBMITTING LJMET -- JERdown"

python -u condor_submitargs.py JERdown

echo "DONE"