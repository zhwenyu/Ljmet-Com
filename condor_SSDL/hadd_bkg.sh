#!/bin/bash

#tt+X
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_TTW.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/TTW/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_TTZ.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/TTZ/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_TTH.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/TTH/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_TTTT.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/TTTT/ | grep '.root'`
#diboson
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_WZ.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/WZ/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_ZZ.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ZZ/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_WpWp.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/WpWp/ | grep '.root'`
##hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_VH.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/VH/ | grep '.root'`
#triboson
#hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_WWW.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/WWW/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_WWZ.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/WWZ/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_WZZ.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/WZZ/ | grep '.root'`
hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_ZZZ.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ZZZ/ | grep '.root'`

#WJets
#hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_WJets_SingleLep.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/WJets-SingleLepSelection/ | grep '.root'`

#qcd
#hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_QCD_Pt15to30_SingleLep.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/QCD-Pt15to30-SingleLepSelection/ | grep '.root'`
#hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_QCD_Pt30to50_SingleLep.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/QCD-Pt30to50-SingleLepSelection/ | grep '.root'`
#hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_QCD_Pt50to80_SingleLep.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/QCD-Pt50to80-SingleLepSelection/ | grep '.root'`
#hadd root://cmseos.fnal.gov//store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/ljmet_trees/ljmet_QCD_Pt80to120_SingleLep.root `xrdfs root://cmseos.fnal.gov ls -u /eos/uscms/store/user/lpctlbsm/clint/Moriond17/25ns/Feb01/QCD-Pt80to120-SingleLepSelection/ | grep '.root'`
