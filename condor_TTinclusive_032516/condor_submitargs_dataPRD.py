import os,sys,datetime

sampleList=[
	'SingleElectron_RRC.txt',
	'SingleMuon_RRC.txt',
	'SingleElectron_PRD_xrd.txt',
	'SingleMuon_PRD_xrd.txt',
]

shift = sys.argv[1]

relBase = '/uscms_data/d3/jmanagan/LJMetSubmit7414'
thisDir = relBase+'/src/LJMet/Com/condor_1lep_120715/' 
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
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt --fileList '+thisDir+'fileLists/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/LJMet_1lepTT_022916 --shift '+shift)
