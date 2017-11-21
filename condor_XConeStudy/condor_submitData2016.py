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
# 	'DoubleMuon_RRH_part1.txt',
# 	'DoubleMuon_RRH_part2.txt',
# 	'DoubleMuon_RRH_part2_Small.txt',

# 	'DoubleEG_RRB.txt',
# 	'DoubleEG_RRC.txt',
# 	'DoubleEG_RRD.txt',
# 	'DoubleEG_RRE.txt',
# 	'DoubleEG_RRF.txt',
# 	'DoubleEG_RRG.txt',
# 	'DoubleEG_RRH.txt',

# 	'MuonEG_RRB.txt',
# 	'MuonEG_RRC.txt',
# 	'MuonEG_RRD.txt',
# 	'MuonEG_RRE.txt',
# 	'MuonEG_RRF.txt',
# 	'MuonEG_RRG.txt',
# 	'MuonEG_RRH.txt',

	'SingleMuon_RRB.txt',
	'SingleMuon_RRC.txt',
	'SingleMuon_RRD.txt',
	'SingleMuon_RRE.txt',
	'SingleMuon_RRF.txt',
	'SingleMuon_RRG.txt',
	'SingleMuon_RRH.txt',

# 	'SingleElectron_RRB.txt',
# 	'SingleElectron_RRC.txt',
# 	'SingleElectron_RRD.txt',
# 	'SingleElectron_RRE.txt',
# 	'SingleElectron_RRF.txt',
# 	'SingleElectron_RRG.txt',
# 	'SingleElectron_RRH.txt',
]

shift = sys.argv[1]

print '====== LJMET SUBMISSION ======'
	
relBase = os.environ['CMSSW_BASE']
print 'Relbase:',relBase

thisDir = relBase+'/src/LJMet/Com/condor_XConeStudy/' 
tarfile = relBase+'.tar'
print 'Making tar:'
if os.path.exists(tarfile):
	print 'tar already exists! Will not re-tar!'
else: 
	os.chdir(relBase)
	os.chdir('../')
	print 'tar --exclude="src/XCone_LJMetSlimmer80x" --exclude="src/LJMet/Com/condor_XConeStudy" --exclude="src/LJMet/Com/condor_TTtrilep" --exclude="src/LJMet/Com/condor_TTinclusive" --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/XCone_LJMetSlimmer80x" --exclude="src/LJMet/Com/condor_XConeStudy" --exclude="src/LJMet/Com/condor_TTtrilep" --exclude="src/LJMet/Com/condor_TTinclusive" --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

cTime=datetime.datetime.now()	
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
# outdir = 'LJMet80x_XCone_2Leps_'+date+'_rizki'
# outdir = 'LJMet80x_XCone_2Leps_2017_3_31_rizki'
# outdir = 'LJMet80x_XCone_1Lep_2017_2017_4_6_rizki'
# outdir = 'LJMet80x_XCone_1Lep_2017_6_15_rizki' #not clustering lep in XCone, properly applied. constituents not saved properly, include area/include ghosts. still use PFchs
# outdir = 'LJMet80x_XCone_2Leps_2017_6_15_rizki' #not clustering lep in XCone, properly applied. constituents not saved properly, include area/include ghosts. still use PFchs
# outdir = 'LJMet80x_XCone_2Leps_2017_6_17_rizki' #not clustering lep in XCone, properly applied.  saving constituents, include area/include ghosts. still use PFchs, include PUPPI, producing XCone based on tauDiff>-30 (PFchs), tauDiff>-25 (PUPPI) 
outdir = 'LJMet80x_XCone_1Lep_2017_6_17_rizki' #not clustering lep in XCone, properly applied.  saving constituents, include area/include ghosts. still use PFchs, include PUPPI, producing XCone based on tauDiff>-30 (PFchs), tauDiff>-25 (PUPPI) 

for sample in sampleList:

	os.system('python condor_submit.py --useMC False --sample '+sample.split('.')[0]+' --json Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --fileList '+thisDir+'fileLists_Mar21/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift)
