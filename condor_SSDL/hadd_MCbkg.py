import os,sys,datetime,time
from ROOT import *
execfile("/uscms_data/d3/rsyarif/EOSSafeUtils.py")

start_time = time.time()

#IO directories must be full paths
input  = sys.argv[1]
output = sys.argv[2]
shift = sys.argv[3]

inputDir='/eos/uscms/store/user/lpcljm/'+input+'/'+shift
outputDir='/eos/uscms/store/user/lpcljm/'+output+'/'+shift

inDir=inputDir[10:]
outDir=outputDir[10:]

os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir)

dirList = [
# 	'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',
# 	'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_combined',
# 	'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',
# 	'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
# 	'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
# 
# # 	'TT_TuneCUETP8M1_13TeV-powheg-pythia8',
# 'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8', #MoriondMC
# # 	'TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8',
# # 	'TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8',
# 'TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8',
# 'TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8',
# 
# 
# #
# # 	'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1',
# # 	'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
# # 	'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
# # 	'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
# # 	'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1',
# #
# 	'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# # 	'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# 	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_combined',
# # 	'DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp',
# 
# 	'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 	'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# 
# 	'DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# 	'DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# 	'DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# 	'DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# 	'DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# 
# #
# 	'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',
# # 	'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# # 	'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# # 	'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# # 	'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# # 	'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# # 	'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# # 	'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# # 	'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
# #
# 	'WW_TuneCUETP8M1_13TeV-pythia8',
# 
# 	'WZ_TuneCUETP8M1_13TeV-pythia8',
# 	'ZZ_TuneCUETP8M1_13TeV-pythia8',
#
# # 	Additionally for multilepton:
# 	'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8',
# 	'ZZTo4L_13TeV_powheg_pythia8',
# 	'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
# 	'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
# 	'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',
# 	'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8',

	#### As requested by Cristina Botta (ARC B2G-17-011).
	'TTGJets',
	'ZGTo2LG',
	'TGJets',
	'WGToLNuG',
	'TTZToLLM1to10',

    ]

# if shift == 'nominal':
# 	dirList.append('DoubleEG_RRB')
# 	dirList.append('DoubleMuon_RRB')
# 	dirList.append('MuonEG_RRB')
# 	dirList.append('DoubleEG_RRC')
# 	dirList.append('DoubleMuon_RRC')
# 	dirList.append('MuonEG_RRC')
# 	dirList.append('DoubleEG_RRD')
# 	dirList.append('DoubleMuon_RRD')
# 	dirList.append('MuonEG_RRD')
# 	dirList.append('DoubleEG_RRE')
# 	dirList.append('DoubleMuon_RRE')
# 	dirList.append('MuonEG_RRE')
# 	dirList.append('DoubleEG_RRF')
# 	dirList.append('DoubleMuon_RRF')
# 	dirList.append('MuonEG_RRF')
# 	dirList.append('DoubleEG_RRG')
# 	dirList.append('DoubleMuon_RRG')
# 	dirList.append('MuonEG_RRG')
# 	dirList.append('DoubleEG_RRH')
# 	dirList.append('DoubleMuon_RRH')
# 	dirList.append('MuonEG_RRH')


for sample in dirList:

    rootfiles = EOSlist_root_files(inputDir+'/'+sample)
#    print 'N root files in',sample,'=',len(rootfiles)
    haddcommand = 'hadd -f root://cmseos.fnal.gov/'+outDir+'/ljmet_'+sample+'_hadd.root '

    print '##########'*15
    print 'HADDING:', sample
    print '##########'*15

    nFilesPerHadd = 1000
#     nFilesPerHadd = 10

    if len(rootfiles) < nFilesPerHadd:
        haddcommand = 'hadd -f root://cmseos.fnal.gov/'+outDir+'/ljmet_'+sample+'_hadd.root '
        for file in rootfiles:
            haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+sample+'/'+file
        os.system(haddcommand)
    else:
        for i in range(int(len(rootfiles)/nFilesPerHadd)):
            haddcommand = 'hadd -f root://cmseos.fnal.gov/'+outDir+'/ljmet_'+sample+'_'+str(i+1)+'_hadd.root '
            begin=i*nFilesPerHadd
            end=begin+nFilesPerHadd
            if len(rootfiles) - end < nFilesPerHadd: end=len(rootfiles)
            for j in range(begin,end):
                haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+sample+'/'+rootfiles[j]
            os.system(haddcommand)
            #        print haddcommand
