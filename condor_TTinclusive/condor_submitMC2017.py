import os,sys,datetime

shift = sys.argv[1]

cernList=[

	# SHOULD BE AT CERN, USE eoscms.cern.ch in condor_submit.py
	# MOVE THINGS HERE TO REMIND YOURSELF WHAT IS AT CERN

	#### TTbar

	]

externalList = [

	# USE cmsxrootd.fnal.gov in condor_submit.py UNLESS YOU CHECK THE LIST LOG FOR CERN PRESENCE	
#	'TprimeTprime_M-1600_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1500_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1300_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1200_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1100_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1000_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1800_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1700_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TprimeTprime_M-1400_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.txt',
#	'TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8.txt',
#	'TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8.txt',
#	'TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8.txt',
	'TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.txt',
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
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/NNKit/.git" --exclude="src/.git" --exclude="src/LJMet-Slimmer" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/NNKit/.git" --exclude="src/.git" --exclude="src/LJMet-Slimmer" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in cernList:
	GenHT = False
	accessor = 'eoscms.cern.ch'
	if 'HT' in sample and 'madgraphMLM' in sample: GenHT = True
	os.system('python condor_submitDeepAK8.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsSummer18/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2018/LJMet94X_1lepTT_081518 --shift '+shift+' --saveGenHT '+str(GenHT)+' --accessor '+accessor)

for sample in externalList:
	GenHT = False
	accessor = 'cmsxrootd.fnal.gov'
	if 'HT' in sample and 'madgraphMLM' in sample: GenHT = True
        os.system('python condor_submitDeepAK8.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsSummer18/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2018/LJMet94X_1lepTT_081518 --shift '+shift+' --saveGenHT '+str(GenHT)+' --accessor '+accessor)

## shift should be (one at a time): nominal, JECup, JECdown, JERup, JERdown
## If you want to use different directory names, edit lines 144 - 147 in condor_submit.py so the config is edited correctly
