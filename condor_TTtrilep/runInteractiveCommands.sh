#!/bin/bash



# ljmet testData_RRB.py >& TriggerDump_RRB.txt
# ljmet testData_RRC.py >& TriggerDump_RRC.txt
# ljmet testData_RRD.py >& TriggerDump_RRD.txt
# ljmet testData_RRE.py >& TriggerDump_RRE.txt
# ljmet testData_RRF.py >& TriggerDump_RRF.txt
# ljmet testData_RRG.py >& TriggerDump_RRG.txt
# ljmet testData_PRH.py >& TriggerDump_PRH.txt

# echo "Checking HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v :"
# grep HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v Trig*.txt
# echo "Checking HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v :"
# grep HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v Trig*.txt
# echo "Checking HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v :"
# grep HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v Trig*.txt
# echo "Checking HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v :"
# grep HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v  Trig*.txt

eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/JECup
eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/JECdown
eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/JERup
eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/JERdown

# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TT_TuneCUETP8M1_13TeV-powheg-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WW_TuneCUETP8M1_13TeV-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/WZ_TuneCUETP8M1_13TeV-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/ZZTo4L_13TeV_powheg_pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/ZZ_TuneCUETP8M1_13TeV-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
# eos root://cmseos.fnal.gov/ rm -r /store/user/lpcljm/LJMet80x_3lepTT_Full2016_2016_12_19_rizki/nominal/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
# 
