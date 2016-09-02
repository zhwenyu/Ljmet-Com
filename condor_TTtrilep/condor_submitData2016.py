import os,sys,datetime

sampleList=[
	# Should be at CERN, use eoscms.cern.ch in condor_submit.py
# 	'SingleElectron_PRB.txt',
# 	'SingleMuon_PRB.txt',
# 	'SingleElectron_PRC.txt',
# 	'SingleMuon_PRC.txt',
# 	'SingleElectron_PRD.txt',
# 	'SingleMuon_PRD.txt',

	'DoubleEG_PRB.txt',
# 	'DoubleMuon_PRB.txt',
# 	'MuonEG_PRB.txt',
	'DoubleEG_PRC.txt',
# 	'DoubleMuon_PRC.txt',
# 	'MuonEG_PRC.txt',
	'DoubleEG_PRD.txt',
# 	'DoubleMuon_PRD.txt',
# 	'MuonEG_PRD.txt',

]

shift = sys.argv[1]

print '====== LJMET SUBMISSION ======'
	
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
outdir = 'LJMet80x_3lepTT_'+date+'_rizki'
# outdir = 'LJMet80x_3lepTT_2016_8_24_rizki'

for sample in sampleList:
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt --fileList '+thisDir+'fileListsSpring16/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift)
