import os,sys,datetime

shift = sys.argv[1]

sampleList=[

	# NOT AT CERN, USE cmsxrootd.fnal.gov in condor_submit.py
	#'TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'X53X53_M-1000_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1200_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1300_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-900_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-180_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-200_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-220_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-250_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-300_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-350_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-400_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-450_13TeV_amcatnlo_pythia8.txt',
# 	'ChargedHiggs_HplusTB_HplusToTB_M-500_13TeV_amcatnlo_pythia8.txt',

	#use xrootd-cms.infn.it/cmsxrootd.fnal.gov:
# 	'TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
# 	'TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 

	# SHOULD BE AT CERN, USE eoscms.cern.ch in condor_submit.py
	'TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  
	'TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  								       
	#'BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	#'BprimeBprime_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  
	#'BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  
	#'BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  	
	#
	#'X53X53_M-1000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1100_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1100_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1300_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1500_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1500_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-1600_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-700_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-700_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	#'X53X53_M-900_LH_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',
	
	]
	
print '=#===== LJMET SUBMISSION ======'
	
relBase = os.environ['CMSSW_BASE']
print 'Relbase:',relBase

thisDir = relBase+'/src/LJMet/Com/condor_TTtrilep/' 
tarfile = relBase+'.tar'
print 'Making tar:'
if os.path.exists(tarfile):
	print 'tar already exists! Will not re-tar!'
else: 
	os.chdir(relBase)
	os.chdir('../')
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
# outdir = 'LJMet80x_3lepTT_'+date+'_rizki'
outdir = 'LJMet80x_3lepTT_2016_8_31_rizki'

for sample in sampleList:
	os.system('python condor_submit.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsSpring16/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift)

## shift should be (one at a time): nominal, JECup, JECdown, JERup, JERdown
## If you want to use different directory names, edit lines 144 - 147 in condor_submit.py so the config is edited correctly
