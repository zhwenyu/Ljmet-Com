import os,sys,datetime

shift = sys.argv[1]

sampleList=[

	# USE cmsxrootd.fnal.gov in condor_submit.py UNLESS YOU CHECK THE LIST LOG FOR CERN PRESENCE
	# SOME OF THESE ARE PROBABLY SPELLED WRONG
	
	# SHOULD BE AT CERN, USE eoscms.cern.ch in condor_submit.py
	# MOVE THINGS HERE TO REMIND YOURSELF WHAT IS AT CERN

	#### Signals
	# 'ChargedHiggs_HplusTB_HplusToTB_M-180_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-200_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-220_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-250_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-3000_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-300_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-350_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-400_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-500_13TeV_amcatnlo_pythia8.txt',
	# 'ChargedHiggs_HplusTB_HplusToTB_M-800_13TeV_amcatnlo_pythia8.txt',

	'BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	
	'TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',

 	'X53X53_M-1000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1000_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1100_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1100_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1200_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1300_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1300_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1500_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1500_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-1600_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-700_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-700_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-900_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	'X53X53_M-900_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	
	#### Drell-yan
	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	
	#### QCD
	'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	
	#### ttbar
	'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',
	'TT_Mtt-700To1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',
	'TT_Mtt-1000ToInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',

	#### single top
	'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.txt',
	'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.txt',
	'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.txt',
	'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',
	'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',

	#### W+Jets
 	'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	# 'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',
	'WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',
	'WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',
	'WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',
	
	# #### Di/Tri boson for CHiggs and multi-lepton
	'WW_TuneCUETP8M1_13TeV-pythia8.txt',
	'WZ_TuneCUETP8M1_13TeV-pythia8.txt',
	'ZZ_TuneCUETP8M1_13TeV-pythia8.txt',

	# 'WWTo2L2Nu_13TeV-powheg.txt',
	# 'WWToLNuQQ_13TeV-powheg.txt',
	
	# 'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt',
	# 'WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8.txt',
	# 'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt',
	# 'WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',
	# 'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.txt',
	# 'WZToLNu2Q_13TeV_powheg_pythia8.txt',
	
	# 'ZZTo2L2Nu_13TeV_powheg_pythia8.txt',
	# 'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt',
	# 'ZZTo2L2Q_13TeV_powheg_pythia8.txt',
	# 'ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8.txt',
	# 'ZZTo2Q2Nu_13TeV_powheg_pythia8.txt',
	# 'ZZTo4L_13TeV-amcatnloFXFX-pythia8.txt',
	# 'ZZTo4L_13TeV_powheg_pythia8.txt',
	
	# 'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	# 'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	# 'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	# 'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	
	# #### Rare SM
	# 'TTTT_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8.txt',
	# 'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.txt',
	# 'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.txt',
	# 'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	# 'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	
	]
	
print '====== LJMET SUBMISSION ======'
	
relBase = os.environ['CMSSW_BASE']
print 'Relbase:',relBase

thisDir = relBase+'/src/LJMet/Com/condor_TTinclusive/' 
tarfile = relBase+'.tar'
print 'Making tar:'
if os.path.exists(tarfile):
	print 'tar already exists! Will not re-tar!'
else: 
	os.chdir(relBase)
	os.chdir('../')
	# YOU NEED TO EXCLUDE ANYTHING ELSE THAT MIGHT LIVE IN THE SAME CMSSW RELEASE
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/PhysicsTools" --exclude="src/EgammaAnalysis" --exclude="src/TopQuarkAnalysis" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/PhysicsTools" --exclude="src/EgammaAnalysis" --exclude="src/TopQuarkAnalysis" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in sampleList:
	GenHT = False
	if 'HT' in sample and 'madgraphMLM' in sample: GenHT = True
	os.system('python condor_submit.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2016/LJMet80X_1lep_022317 --shift '+shift+' --saveGenHT '+str(GenHT))

## shift should be (one at a time): nominal, JECup, JECdown, JERup, JERdown
## If you want to use different directory names, edit lines 144 - 147 in condor_submit.py so the config is edited correctly
