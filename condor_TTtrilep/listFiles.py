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
   '/TprimeTprime_M-1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM',


# 	#WZTo3LNu
# 	'/WZTo3LNu_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
# 	#ZZTo4L
# 	'/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',
# 	'/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM',
# 	#WWW 4F
# 	'/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
# 	#WWZ
# 	'/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
# 	#WZZ
# 	'/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',
# 	#ZZZ
# 	'/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
# 	#TTWJetsToLNu
# 	'/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',
# 	#TTZToLLNuNu M-10
# 	'/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',
# 	
#	#WW
#	'/WW_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',

	#DYJetsToLL M-50 --> Julie has this	
	#WJetsToLNu --> Julie has this
	#TT --> Julis has this

    ]
# 
# datalistRRB = [
# '/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD',
# '/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD',
# '/MuonEG/Run2017B-31Mar2018-v1/MINIAOD',
# ]
# datalistRRC = [
# '/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD',
# '/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD',
# '/MuonEG/Run2017C-31Mar2018-v1/MINIAOD',
# ]
# datalistRRD = [
# '/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD',
# '/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD',
# '/MuonEG/Run2017D-31Mar2018-v1/MINIAOD',
# ]
# datalistRRE = [
# '/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD',
# '/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD',
# '/MuonEG/Run2017E-31Mar2018-v1/MINIAOD',
# ]
# datalistRRF = [
# '/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD',
# '/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD',
# '/MuonEG/Run2017F-31Mar2018-v1/MINIAOD',
# ]
# datalistRRF_v2 = [
# '/DoubleMuon/Run2017F-09May2018-v1/MINIAOD',
# '/DoubleEG/Run2017F-09May2018-v1/MINIAOD',
# '/MuonEG/Run2017F-09May2018-v1/MINIAOD',
# ]
# 

for sample in samplelist:
    print 'listing files in',sample
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" &>> fileLists_Nov2-2018/'+sample.split('/')[1]+'.txt')
    #os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')
    #os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="dataset = '+sample+' | grep dataset.nevents" --key ~/.globus/userkey.pem --cert ~/.globus/usercert.pem')

# for sample in datalistRRB:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_Nov2-2018/'+sample.split('/')[1]+'_RRB.txt')
# 
# for sample in datalistRRC:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_Nov2-2018/'+sample.split('/')[1]+'_RRC.txt')
# 
# for sample in datalistRRD:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_Nov2-2018/'+sample.split('/')[1]+'_RRD.txt')
# 
# for sample in datalistRRE:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_Nov2-2018/'+sample.split('/')[1]+'_RRE.txt')
# 
# for sample in datalistRRF:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_Nov2-2018/'+sample.split('/')[1]+'_RRF.txt')
# 
# for sample in datalistRRF_v2:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileLists_Nov2-2018/'+sample.split('/')[1]+'_RRF_v2.txt')
# 
# 
