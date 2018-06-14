import os,sys

samplelist = [

   '/TprimeTprime_M-1800_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1600_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1400_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1100_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM',
   '/TprimeTprime_M-1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'

    ]

datalistRRC = [
#    '/SingleElectron/Run2015C_25ns-05Oct2015-v1/MINIAOD',
#    '/SingleMuon/Run2015C_25ns-05Oct2015-v1/MINIAOD',
#    '/DoubleEG/Run2015C_25ns-05Oct2015-v1/MINIAOD',
#    '/DoubleMuon/Run2015C_25ns-05Oct2015-v1/MINIAOD',
#    '/MuonEG/Run2015C_25ns-05Oct2015-v1/MINIAOD',
]

datalistRRD = [
#    '/SingleElectron/Run2015D-05Oct2015-v1/MINIAOD',
#    '/SingleMuon/Run2015D-05Oct2015-v1/MINIAOD',
#    '/DoubleMuon/Run2015D-05Oct2015-v1/MINIAOD',
#    '/MuonEG/Run2015D-05Oct2015-v1/MINIAOD',
#    '/DoubleEG/Run2015D-05Oct2015-v1/MINIAOD',
]

datalistPRD = [
#    '/SingleElectron/Run2015D-PromptReco-v4/MINIAOD',
#    '/SingleMuon/Run2015D-PromptReco-v4/MINIAOD',
#    '/DoubleMuon/Run2015D-PromptReco-v4/MINIAOD',
#    '/MuonEG/Run2015D-PromptReco-v4/MINIAOD',
#    '/DoubleEG/Run2015D-PromptReco-v4/MINIAOD',
]

datalistPRE = [
    '/SingleMuon/Run2016E-PromptReco-v2/MINIAOD',
    '/SingleElectron/Run2016E-PromptReco-v2/MINIAOD',
]

datalistPRF = [
    '/SingleMuon/Run2016F-PromptReco-v1/MINIAOD',
    '/SingleElectron/Run2016F-PromptReco-v1/MINIAOD',
]

datalistPRG = [
    '/SingleMuon/Run2016G-PromptReco-v1/MINIAOD',
    '/SingleElectron/Run2016G-PromptReco-v1/MINIAOD',
]

datalistRR = [
    '/JetHT/Run2016B-23Sep2016-v3/MINIAOD',
    '/JetHT/Run2016C-23Sep2016-v1/MINIAOD',
    '/JetHT/Run2016D-23Sep2016-v1/MINIAOD',
    '/JetHT/Run2016E-23Sep2016-v1/MINIAOD',
    '/JetHT/Run2016F-23Sep2016-v1/MINIAOD',
    '/JetHT/Run2016G-23Sep2016-v1/MINIAOD',
    ]

datalistPRH = [
    '/JetHT/Run2016H-PromptReco-v2/MINIAOD',
    '/JetHT/Run2016H-PromptReco-v3/MINIAOD',
]

for sample in samplelist:
    print 'listing files in',sample
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileListsSummer18/'+sample.split('/')[1]+'.txt')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')
    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="dataset = '+sample+' | grep dataset.nevents" ')

# for sample in datalistPRE:
#    print 'listing files in',sample
#    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileListsSpring16/'+sample.split('/')[1]+'_PRE.txt')
#    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')

# for sample in datalistPRF:
#    print 'listing files in',sample
#    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileListsSpring16/'+sample.split('/')[1]+'_PRF.txt')
#    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')

# for sample in datalistPRG:
#    print 'listing files in',sample
#    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >& fileListsSpring16/'+sample.split('/')[1]+'_PRG.txt')
#    os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')

# for sample in datalistPRH:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >> fileListsSpring16/'+sample.split('/')[1]+'_PRH.txt')
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')

# for sample in datalistRR:
#     print 'listing files in',sample
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="file dataset = '+sample+'" >> fileListsSpring16/'+sample.split('/')[1]+'_RRBCDEFG.txt')
#     os.system('/cvmfs/cms.cern.ch/common/dasgoclient --limit=0 --query="site dataset = '+sample+'" | grep "T2_CH_CERN" ')

