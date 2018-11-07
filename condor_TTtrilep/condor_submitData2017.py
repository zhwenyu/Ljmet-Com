import os,sys,datetime

sampleList=[
	# Should be at CERN, use eoscms.cern.ch in condor_submit.py

	'DoubleMuon_RRB.txt',
	'DoubleMuon_RRC.txt',
	'DoubleMuon_RRD.txt',
	'DoubleMuon_RRE.txt',
	'DoubleMuon_RRF.txt',
	'DoubleMuon_RRF_v2.txt',

	'DoubleEG_RRB.txt',
	'DoubleEG_RRC.txt',
	'DoubleEG_RRD.txt',
	'DoubleEG_RRE.txt',
	'DoubleEG_RRF.txt',
	'DoubleEG_RRF_v2.txt',

	'MuonEG_RRB.txt',
	'MuonEG_RRC.txt',
	'MuonEG_RRD.txt',
	'MuonEG_RRE.txt',
	'MuonEG_RRF.txt',
	'MuonEG_RRF_v2.txt',

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
# 	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" --exclude="*.scram" -zcf '+tarfile+' '+relBase.split('/')[-1]+'/'
# 	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" --exclude="*.scram" -zcf '+tarfile+' '+relBase.split('/')[-1])
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" --exclude="*.scram" -zcf '+tarfile+' '+relBase.split('/')[-1]+'/src/*' #only tar on the level where LJMet is
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" --exclude="*.scram" -zcf '+tarfile+' '+relBase.split('/')[-1]+'/src/*') #only tar on the level where LJMet is
	os.chdir(thisDir)

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
outdir = 'LJMet94x_3lepTT_2017datasets_'+date+'_rizki'

for sample in sampleList:
	accessor = 'cmsxrootd.fnal.gov'
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt --fileList '+thisDir+'fileLists_Nov2-2018/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift+' --accessor '+accessor)
