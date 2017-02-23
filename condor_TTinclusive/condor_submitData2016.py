import os,sys,datetime

sampleList=[
	# Should be at CERN, use eoscms.cern.ch in condor_submit.py
	#'JetHT_PRH.txt',
	#'JetHT_RRBCDEFG.txt',
	'SingleElectron_RRBCDEFGH.txt',
	'SingleMuon_RRBCDEFGH.txt',

	# At LPC, use cmseos.fnal.gov in condor_submit.py
	#'SingleElectron_PRH.txt',
	#'SingleMuon_PRH.txt',
	#'SingleElectron_RRBCDEFG.txt',
	#'SingleMuon_RRBCDEFG.txt',
]

shift = sys.argv[1]

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
	# YOU NEED TO EXCLUDE EVERYTHING THAT MIGHT BE IN THE SAME CMSSW RELEASE
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/PhysicsTools" --exclude="src/EgammaAnalysis" --exclude="src/TopQuarkAnalysis" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/PhysicsTools" --exclude="src/EgammaAnalysis" --exclude="src/TopQuarkAnalysis" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in sampleList:
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2016/LJMet80X_1lep_022317 --shift '+shift)
