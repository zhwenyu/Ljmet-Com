import os,sys,datetime

sampleList=[
	# Should be at CERN, use eoscms.cern.ch in condor_submit.py
	#'JetHT_PRH.txt',
	#'JetHT_RRBCDEFG.txt',
	#'SingleElectron_RRBCDEFGH.txt',
	'SingleElectron_RRHv3.txt',
	'SingleMuon_RRHv3.txt',
	#'SingleMuon_RRBCDEF.txt',
	#'DoubleEG_PRRBCDEFGH.txt',	
	#'DoubleMuon_PRRBCDEFGH.txt',	
	#'MuonEG_PRRBCDEFGH.txt',	
	#'JetHT_RRCDEFG.txt',

	# At LPC, use cmseos.fnal.gov in condor_submit.py
	#'SingleElectron_PRH.txt',
	#'SingleMuon_PRH.txt',
	#'SingleElectron_RRBCDEFG.txt',
	#'SingleMuon_RRBCDEFG.txt',
]

notcernList=[
	#'SingleMuon_RRGH.txt',
	#'JetHT_RRBH.txt',
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
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/PhysicsTools" --exclude="src/EgammaAnalysis" --exclude="src/TopQuarkAnalysis" --exclude="src/LJMetSlimmer" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/PhysicsTools" --exclude="src/EgammaAnalysis" --exclude="src/TopQuarkAnalysis" --exclude="src/LJMetSlimmer" --exclude="src/X53ThirteenTeVAnalysisCode" --exclude="src/X53ThirteenTeVAnalysisCode.tar" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in sampleList:
	accessor = 'eoscms.cern.ch'
	#os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2016/LJMet80X_1lep_041917 --shift '+shift+' --accessor '+accessor)
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2016/LJMet80X_1lep_041917 --shift '+shift+' --accessor '+accessor)

for sample in notcernList:
	accessor = 'cmsxrootd.fnal.gov'
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2016/LJMet80X_1lep_041917 --shift '+shift+' --accessor '+accessor)
