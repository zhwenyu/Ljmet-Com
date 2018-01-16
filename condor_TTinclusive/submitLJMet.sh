#!/bin/bash

echo "SUBMITTING LJMET -- PR 2016R,C,D Data"

#cp TTsinglelep_DataElCheck_cfg.py ljmet_cfg.py
cp TTsinglelep_Data_cfg.py ljmet_cfg.py

python -u condor_submitData2016.py nominal #>& submit_Data.log &

echo "SUBMITTING LJMET -- nominal"

cp TTsinglelep_MC_cfg.py ljmet_cfg.py

python -u condor_submitMoriond2017.py nominal #>& submit_MC_nominal.log &

echo "SUBMITTING LJMET -- JECup"

python -u condor_submitMoriond2017.py JECup #>& submit_MC_JECup.log &

echo "SUBMITTING LJMET -- JECdown"

python -u condor_submitMoriond2017.py JECdown #>& submit_MC_JECdown.log &

echo "SUBMITTING LJMET -- JERup"

python -u condor_submitMoriond2017.py JERup #>& submit_MC_JERup.log &

echo "SUBMITTING LJMET -- JERdown"

python -u condor_submitMoriond2017.py JERdown #>& submit_MC_JERdown.log &

echo "DONE"
