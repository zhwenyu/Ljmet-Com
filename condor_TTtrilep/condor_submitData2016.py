import os,sys,datetime

sampleList=[
	# Should be at CERN, use eoscms.cern.ch in condor_submit.py

# 	'DoubleMuon_RRB.txt',
# 	'DoubleMuon_RRC.txt',
# 	'DoubleMuon_RRD.txt',
# 	'DoubleMuon_RRE.txt',
# 	'DoubleMuon_RRF.txt',
# 	'DoubleMuon_RRG.txt',
# 	'DoubleMuon_RRH.txt',

# 	'DoubleEG_RRB.txt',
# 	'DoubleEG_RRC.txt',
# 	'DoubleEG_RRD.txt',
# 	'DoubleEG_RRE.txt',
# 	'DoubleEG_RRF.txt',
# 	'DoubleEG_RRG.txt',
	'DoubleEG_RRH.txt',

# 	'MuonEG_RRB.txt',
# 	'MuonEG_RRC.txt',
# 	'MuonEG_RRD.txt',
# 	'MuonEG_RRE.txt',
# 	'MuonEG_RRF.txt',
# 	'MuonEG_RRG.txt',
# 	'MuonEG_RRH.txt',

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
# outdir = 'LJMet80x_3lepTT_Full2016_Moriond17_reMiniAOD_nuBTVSF_modMETfilt_'+date+'_rizki'
# outdir = 'LJMet80x_3lepTT_Full2016_Moriond17_reMiniAOD_nuBTVSF_modMETfilt_newRunH_'+date+'_rizki'
outdir = 'LJMet80x_3lepTT_Full2016_Moriond17_reMiniAOD_nuBTVSF_modMETfilt_newRunH_2017_3_21_rizki'

for sample in sampleList:
# 	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --fileList '+thisDir+'fileLists_Feb24/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift)
	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --fileList '+thisDir+'fileLists_Mar21/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift)
