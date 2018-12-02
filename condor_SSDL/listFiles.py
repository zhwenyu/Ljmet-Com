import os,sys

samplelist = [

	# TT
#    '/TprimeTprime_M-1800_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1600_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1400_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1100_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
#    '/TprimeTprime_M-1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM',


# 	#WZTo3LNu
	'/WZTo3LNu_TuneCP5_13TeV-powheg-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	#ZZTo4L
	'/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	#WWW 4F
	'/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	#WWZ
	'/WWZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	#WZZ
	'/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	#ZZZ
	'/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	#TTWJetsToLNu
	'/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	#TTZToLLNuNu M-10
	'/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',
# 	
#	#WW
#	'/WW_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',

	#DYJetsToLL M-50 --> Julie has this	
	#WJetsToLNu --> Julie has this
	#TT --> Julis has this
	
	#TTTT
	'/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12-v1/MINIAODSIM',

    ]

datalistRunA = [
# '/EGamma/Run2018A-22Jun2018-v1/MINIAOD',
# '/DoubleMuon/Run2018A-17Sep2018-v2/MINIAOD',
# '/MuonEG/Run2018A-PromptReco-v3/MINIAOD'
]
datalistRunB = [
# '/EGamma/Run2018B-17Sep2018-v1/MINIAOD',
# '/EGamma/Run2018B-26Sep2018-v1/MINIAOD',
# '/DoubleMuon/Run2018B-17Sep2018-v1/MINIAOD',
# '/MuonEG/Run2018B-PromptReco-v1/MINIAOD', #chose this because it had more files
]
datalistRunC = [
# '/EGamma/Run2018C-17Sep2018-v1/MINIAOD',
# '/DoubleMuon/Run2018C-17Sep2018-v1/MINIAOD',
# '/MuonEG/Run2018C-PromptReco-v3/MINIAOD',
]
datalistRunD = [
# '/EGamma/Run2018D-PromptReco-v2/MINIAOD',
# '/DoubleMuon/Run2018D-PromptReco-v2/MINIAOD',
# '/MuonEG/Run2018D-PromptReco-v2/MINIAOD',
]


for sample in samplelist:
    print 'listing files in',sample
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" &>> fileLists_2018dataset_Nov21-2018/'+sample.split('/')[1]+'.txt')
    #os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')
    #os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="dataset = '+sample+' | grep dataset.nevents" --key ~/.globus/userkey.pem --cert ~/.globus/usercert.pem')

for sample in datalistRunA:
    print 'listing files in',sample
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_2018dataset_Nov21-2018/'+sample.split('/')[1]+'_RunA.txt')

for sample in datalistRunB:
    print 'listing files in',sample
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_2018dataset_Nov21-2018/'+sample.split('/')[1]+'_RunB.txt')

for sample in datalistRunC:
    print 'listing files in',sample
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_2018dataset_Nov21-2018/'+sample.split('/')[1]+'_RunC.txt')

for sample in datalistRunD:
    print 'listing files in',sample
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_2018dataset_Nov21-2018/'+sample.split('/')[1]+'_RunD.txt')

