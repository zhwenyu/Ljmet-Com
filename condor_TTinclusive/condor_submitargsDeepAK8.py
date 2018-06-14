import os,sys,datetime

shift = sys.argv[1]

cernList=[

	# SHOULD BE AT CERN, USE eoscms.cern.ch in condor_submit.py
	# MOVE THINGS HERE TO REMIND YOURSELF WHAT IS AT CERN

	#### TTbar
	#'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.txt',
	#'TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.txt',
	'TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',
	'TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',

	# #### ISR/FSR
# 	'TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8.txt',
# 	'TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8.txt',
# 	'TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8.txt',
# 	'TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8.txt',
	
	]

externalList = [

	# USE cmsxrootd.fnal.gov in condor_submit.py UNLESS YOU CHECK THE LIST LOG FOR CERN PRESENCE	
	# 'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',
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
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in cernList:
	GenHT = False
	accessor = 'eoscms.cern.ch'
	if 'HT' in sample and 'madgraphMLM' in sample: GenHT = True
	os.system('python condor_submitDeepAK8.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/CHiggs/LJMet80X_1lepTTsplit_041117 --shift '+shift+' --saveGenHT '+str(GenHT)+' --accessor '+accessor)

for sample in externalList:
	GenHT = False
	accessor = 'cmsxrootd.fnal.gov'
	if 'HT' in sample and 'madgraphMLM' in sample: GenHT = True
	os.system('python condor_submitDeepAK8.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/CHiggs/LJMet80X_1lepTTsplit_041117 --shift '+shift+' --saveGenHT '+str(GenHT)+' --accessor '+accessor)

## shift should be (one at a time): nominal, JECup, JECdown, JERup, JERdown
## If you want to use different directory names, edit lines 144 - 147 in condor_submit.py so the config is edited correctly
