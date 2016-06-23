#!/bin/tcsh -f

#foreach mass ( 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 )
#   echo "Dumping /WprimeToTB_TToLep_M-${mass}_RH_TuneCUETP8M1_13TeV-comphep-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM:"
#   set vers = 1
#   if (mass == 4000 || mass == 3500) then
#     set vers = 3
#   endif
#   eval das_client --query=\"file dataset=/WprimeToTB_TToLep_M-${mass}_RH_TuneCUETP8M1_13TeV-comphep-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM\" --format=plain --limit=0 > python/Samples_2016/WprimeToTB_TToLep_M_${mass}_RH_TuneCUETP8M1_13TeV_comphep_pythia8_RunIISpring16MiniAODv2_PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0_v${vers}_cff.txt 
#   echo 'Done!'
#end

#foreach dataset ( TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext3-v1 WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1 DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1 WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1 WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1 QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1 )
#   set tmp1 = `echo ${dataset} | sed 's/-/_/g'`
#   set tmp2 = `echo ${tmp1} | sed 's/\//_/g'`
#   echo "Dumping /${dataset}/MINIAODSIM:"
#   eval das_client --query=\"file dataset=/${dataset}/MINIAODSIM\" --limit=0 > python/Samples_2016/${tmp2}_cff.txt 
#   echo 'Done!'
#end

foreach dataset ( Run2016B-PromptReco-v2 )
   set tmp = `echo ${dataset} | sed 's/-/_/g'`
   echo "Dumping /SingleLepton/${dataset}/MINIAOD:"
   eval das_client --query=\"file dataset=/SingleElectron/${dataset}/MINIAOD\" --limit=0 > python/Samples_2016/SingleElectron_${tmp}_cff.txt 
   eval das_client --query=\"file dataset=/SingleMuon/${dataset}/MINIAOD\" --limit=0 > python/Samples_2016/SingleMuon_${tmp}_cff.txt 
   echo 'Done!'
end
