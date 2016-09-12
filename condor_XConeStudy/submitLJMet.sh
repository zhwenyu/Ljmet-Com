#!/bin/bash

echo "SUBMITTING LJMET -- PR 2016R,C,D Data"

cp TTtrilep_Data_cfg.py ljmet_cfg.py

# python -u condor_submitData2016.py nominal >& submit_Data.log &
python -u condor_submitData2016.py nominal 

# echo "SUBMITTING LJMET -- nominal MC"
# 
# cp TTtrilep_MC_cfg.py ljmet_cfg.py
# 
# # python -u condor_submitSpring2016.py nominal >& submit_MC_nominal.log &
# python -u condor_submitSpring2016.py nominal

#echo "SUBMITTING LJMET -- JECup"

#python -u condor_submitSpring2016.py JECup #>& submit_MC_JECup.log &

#echo "SUBMITTING LJMET -- JECdown"

#python -u condor_submitSpring2016.py JECdown #>& submit_MC_JECdown.log &

#echo "SUBMITTING LJMET -- JERup"

#python -u condor_submitSpring2016.py JERup #>& submit_MC_JERup.log &

#echo "SUBMITTING LJMET -- JERdown"

#python -u condor_submitSpring2016.py JERdown #>& submit_MC_JERdown.log &

#echo "SUBMITTING LJMET -- nominal"

#cp TTtrilep_reHLT_cfg.py ljmet_cfg.py

#python -u condor_submitReHLT.py nominal #>& submit_MC_nominal.log &

#echo "SUBMITTING LJMET -- JECup"

#python -u condor_submitReHLT.py JECup #>& submit_MC_JECup.log &

#echo "SUBMITTING LJMET -- JECdown"

#python -u condor_submitReHLT.py JECdown #>& submit_MC_JECdown.log &

#echo "SUBMITTING LJMET -- JERup"

#python -u condor_submitReHLT.py JERup #>& submit_MC_JERup.log &

#echo "SUBMITTING LJMET -- JERdown"

#python -u condor_submitReHLT.py JERdown #>& submit_MC_JERdown.log &

echo "DONE"
