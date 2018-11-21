#!/bin/bash

rm -v /uscms_data/d3/rsyarif/Brown2018/TT_BB_trilepton_LJMET_2017data/CMSSW_9_4_11.tar 


echo "SUBMITTING LJMET -- RR 2017Data"

cp TTtrilep_Data_cfg.py ljmet_cfg.py

python -u condor_submitData2017.py nominal | tee submit_Data2017.log 

# echo "SUBMITTING LJMET -- nominal MC"
# 
# cp -v TTtrilep_MC_cfg.py ljmet_cfg.py
# 
# python -u condor_submitMC2017.py nominal | tee submit_MC2017_nominal.log 

# echo "SUBMITTING LJMET -- JECup"
# 
# python -u condor_submitMC2017.py JECup #>& submit_MC2017_JECup.log &
# 
# echo "SUBMITTING LJMET -- JECdown"
# 
# python -u condor_submitMC2017.py JECdown #>& submit_MC2017_JECdown.log &
# 
# echo "SUBMITTING LJMET -- JERup"
# 
# python -u condor_submitMC2017.py JERup #>& submit_MC2017_JERup.log &
# 
# echo "SUBMITTING LJMET -- JERdown"
# 
# python -u condor_submitMC2017.py JERdown #>& submit_MC2017_JERdown.log &


echo "DONE"
