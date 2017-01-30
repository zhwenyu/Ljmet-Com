#!/bin/bash

rm -v /uscms_data/d3/rsyarif/FermilabFall2016/produceLJMetNtuple2016_FullDataset_v2/CMSSW_8_0_20.tar 

cd ..
# source switchToTpTp.sh
scramv1 b -j10
cd -

echo "SUBMITTING LJMET -- RR 2016B-G PR 2016H Data"

cp TTtrilep_Data_cfg.py ljmet_cfg.py

python -u condor_submitData2016.py nominal #>& submit_Data.log &

echo "SUBMITTING LJMET -- nominal MC"

cp -v TTtrilep_MC_cfg.py ljmet_cfg.py

python -u condor_submitMoriond2017.py nominal #>& submit_MC_nominal.log &

echo "SUBMITTING LJMET -- JECup"

python -u condor_submitMoriond2017.py JECup #>& submit_MC_JECup.log &

echo "SUBMITTING LJMET -- JECdown"

python -u condor_submitMoriond2017.py JECdown #>& submit_MC_JECdown.log &

echo "SUBMITTING LJMET -- JERup"

python -u condor_submitMoriond2017.py JERup #>& submit_MC_JERup.log &

echo "SUBMITTING LJMET -- JERdown"

python -u condor_submitMoriond2017.py JERdown #>& submit_MC_JERdown.log &
# 

echo "DONE"
