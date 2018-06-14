#!/bin/bash

echo "SUBMITTING LJMET -- nominal"
cp TTsinglelep_MCDeepAK8_cfg.py ljmet_cfg.py
python -u condor_submitargsDeepAK8.py nominal

echo "SUBMITTING LJMET -- JECup"
cp TTsinglelep_MCDeepAK8_cfg.py ljmet_cfg.py
python -u condor_submitargsDeepAK8.py JECup

echo "SUBMITTING LJMET -- JECdown"
cp TTsinglelep_MCDeepAK8_cfg.py ljmet_cfg.py
python -u condor_submitargsDeepAK8.py JECdown

echo "SUBMITTING LJMET -- JERup"
cp TTsinglelep_MCDeepAK8_cfg.py ljmet_cfg.py
python -u condor_submitargsDeepAK8.py JERup

echo "SUBMITTING LJMET -- JERdown"
cp TTsinglelep_MCDeepAK8_cfg.py ljmet_cfg.py
python -u condor_submitargsDeepAK8.py JERdown

echo "DONE"
