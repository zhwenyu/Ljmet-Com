import os,sys,datetime

sampleList=[
	# Should be at CERN, use eoscms.cern.ch in condor_submit.py

# 	'DoubleMuon_RunA.txt',
# 	'DoubleMuon_RunB.txt',
# 	'DoubleMuon_RunC.txt',
# 	'DoubleMuon_RunD.txt',
# 
# 
# 	'EGamma_RunA.txt',
	'EGamma_RunB_26Sep2018.txt',
# 	'EGamma_RunC.txt',
# 	'EGamma_RunD.txt',

# 	'EGamma_RunB.txt',
# 
# 	'MuonEG_RunA.txt',
# 	'MuonEG_RunB.txt',
# 	'MuonEG_RunC.txt',
# 	'MuonEG_RunD.txt',

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
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" --exclude="*.SCRAM" -zcf '+tarfile+' '+relBase.split('/')[-1]+'/src/*' #only tar on the level where LJMet is
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" --exclude="*.SCRAM" -zcf '+tarfile+' '+relBase.split('/')[-1]+'/src/*') #only tar on the level where LJMet is
	os.chdir(thisDir)

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
# outdir = 'LJMet102x_3lepTT_2018datasets_'+date+'_rizki'
outdir = 'LJMet102x_3lepTT_2018datasets_2018_11_22_rizki'

for sample in sampleList:
	accessor = 'cmsxrootd.fnal.gov'
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt --fileList '+thisDir+'fileLists_2018dataset_Nov21-2018/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift+' --accessor '+accessor)
