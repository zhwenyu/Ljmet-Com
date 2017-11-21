#!/bin/bash

rm -vi /uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/*.tar

###dilep for Zjets
# 
# echo "SUBMITTING LJMET (dilep)-- RR 2016B-H Data"
# 
# cp -v XCone_Data_dilep_cfg.py ljmet_cfg.py
# 
# python -u condor_submitData2016.py nominal #>& submit_Data.log &
# 
# echo "SUBMITTING LJMET (MC dilep)-- nominal MC"
# 
# cp -v XCone_MC_dilep_cfg.py ljmet_cfg.py
# 
# python -u condor_submitMoriond2017.py nominal #>& submit_MC_nominal.log &


###singlelep for QCD/ttbar
# 
# echo "SUBMITTING LJMET (singlelep)-- RR 2016B-H Data"
# 
# cp -v  XCone_Data_singlelep_cfg.py ljmet_cfg.py
# 
# python -u condor_submitData2016.py nominal #>& submit_Data.log &

echo "SUBMITTING LJMET (MC singlelep)-- nominal MC"

cp -v  XCone_MC_singlelep_cfg.py ljmet_cfg.py

python -u condor_submitMoriond2017.py nominal #>& submit_MC_nominal.log &


### no leps - ttbar && TpTp
# 
# echo "SUBMITTING LJMET -- nominal MC"
# 
# cp -v  XCone_MC_cfg.py ljmet_cfg.py
# 
# # python -u condor_submitMoriond2017.py nominal >& submit_MC_nominal.log &
# python -u condor_submitMoriond2017.py nominal

#echo "SUBMITTING LJMET -- JECup"

#python -u condor_submitMoriond2016.py JECup #>& submit_MC_JECup.log &

#echo "SUBMITTING LJMET -- JECdown"

#python -u condor_submitMoriond2016.py JECdown #>& submit_MC_JECdown.log &

#echo "SUBMITTING LJMET -- JERup"

#python -u condor_submitMoriond2016.py JERup #>& submit_MC_JERup.log &

#echo "SUBMITTING LJMET -- JERdown"

#python -u condor_submitMoriond2016.py JERdown #>& submit_MC_JERdown.log &

echo "DONE"
