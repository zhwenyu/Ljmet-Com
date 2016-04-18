import os,sys,datetime

sampleList=[
	'SingleElectron_RRC.txt',
	'SingleMuon_RRC.txt',
	'SingleElectron_PRD_xrd.txt',
	'SingleMuon_PRD_xrd.txt',
]

shift = sys.argv[1]

relBase = os.environ['CMSSW_BASE']
thisDir = relBase+'/src/LJMet/Com/condor/' 
tarfile = relBase+'.tar'
if os.path.exists(tarfile):
	print 'tar already exists! Will not re-tar!'
else: 
	os.chdir(relBase)
	os.chdir('../')
	print 'tar -cf',tarfile,' ',relBase.split('/')[-1],'/'
	os.system('tar -cf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in sampleList:
	os.system('python TTsinglelep_condorsubmit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt --fileList '+thisDir+'fileLists/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/LJMet_1lepTT_020516 --shift '+shift)
