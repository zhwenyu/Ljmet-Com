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

# signalList = [
# #    'TprimeTprime_M-700',
#     'TprimeTprime_M-800',
#     'TprimeTprime_M-900',
#     'TprimeTprime_M-1000',
#     'TprimeTprime_M-1100',
#     'TprimeTprime_M-1200',
#     'TprimeTprime_M-1300',
#     'TprimeTprime_M-1400',
#     'TprimeTprime_M-1500',
#     'TprimeTprime_M-1600',
#     'TprimeTprime_M-1700',
#     'TprimeTprime_M-1800',
#     ]
# 
# signalOutList = ['BWBW','TZBW','THBW','TZTH','TZTZ','THTH']
# 
# for sample in signalList:
# #     for outlabel in signalOutList:
# 
#         rootfiles = EOSlist_root_files(inputDir+'/'+sample+'/')
#         print 'N root files in',sample,'=',len(rootfiles)
#         haddcommand = 'hadd -f root://cmseos.fnal.gov/'+outDir+'/ljmet_'+sample+'_hadd.root '
# 
#         print '##########'*15
#         print 'HADDING:', sample
#         print '##########'*15
# 
#         for file in rootfiles:
#            
#             haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+sample+'/'+file
#         
#         os.system(haddcommand)
#         #print haddcommand

signalList = [
   #'BprimeBprime_M-700',
    'BprimeBprime_M-800',
    'BprimeBprime_M-900',
    'BprimeBprime_M-1000',
    'BprimeBprime_M-1100',
    'BprimeBprime_M-1200',
    'BprimeBprime_M-1300',
    'BprimeBprime_M-1400',
    'BprimeBprime_M-1500',
    'BprimeBprime_M-1600',
    'BprimeBprime_M-1700',
    'BprimeBprime_M-1800',
    ]

signalOutList = ['TWTW','BZTW','BHTW','BZBH','BZBZ','BHBH']

for sample in signalList:
#     for outlabel in signalOutList:

        rootfiles = EOSlist_root_files(inputDir+'/'+sample+'/')
#        print 'N root files in',sample,'=',len(rootfiles)
        haddcommand = 'hadd -f root://cmseos.fnal.gov/'+outDir+'/ljmet_'+sample+'_hadd.root '

        print '##########'*15
        print 'HADDING:', sample
        print '##########'*15

        for file in rootfiles:

            haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+sample+'/'+file

        os.system(haddcommand)
#        print haddcommand

