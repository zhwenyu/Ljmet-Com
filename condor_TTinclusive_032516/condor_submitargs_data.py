import os,sys,datetime

shift = sys.argv[1]

sampleList=[
'DoubleEG_Run2015C_16Dec2015.txt',
'DoubleEG_Run2015D_16Dec2015.txt',
'DoubleMuon_Run2015C_16Dec2015.txt',
'DoubleMuon_Run2015D_16Dec2015.txt',
'MuonEG_Run2015C_16Dec2015.txt',
'MuonEG_Run2015D_16Dec2015.txt',
# 	'SingleElectron_Run2015C_16Dec2015.txt',
# 	'SingleElectron_Run2015D_16Dec2015.txt',
# 	'SingleMuon_Run2015C_16Dec2015.txt',
# 	'SingleMuon_Run2015D_16Dec2015.txt',
]

print '====== LJMET SUBMISSION ======'
	
relBase = os.environ['CMSSW_BASE']
print 'Relbase:',relBase

thisDir = relBase+'/src/LJMet/Com/condor_TTinclusive_032516/' 
tarfile = relBase+'.tar'
print 'Making tar:'
if os.path.exists(tarfile):
	print 'tar already exists! Will not re-tar!'
else: 
	os.chdir(relBase)
	os.chdir('../')
	print 'tar --exclude="src/LJMet/Com/.git" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in sampleList:
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt --fileList '+thisDir+'fileLists/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/LJMet_3lepTT_051316 --shift '+shift)
							
## shift should be (one at a time): nominal, JECup, JECdown, JERup, JERdown
## If you want to use different directory names, edit lines 144 - 147 in condor_submit.py so the config is edited correctly

