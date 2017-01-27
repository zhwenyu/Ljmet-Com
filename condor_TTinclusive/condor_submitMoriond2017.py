import os,sys,datetime

shift = sys.argv[1]

sampleList=[

	# USE cmsxrootd.fnal.gov in condor_submit.py UNLESS YOU CHECK THE LIST LOG FOR CERN PRESENCE
	# SOME OF THESE ARE PROBABLY SPELLED WRONG
	'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',				  
										  
	'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.txt',	  
	'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',	  
	'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',	  
	'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt', 
	'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',	  
	
	'DYJetsToLL_M-50_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'DYJetsToLL_M-50_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'DYJetsToLL_M-50_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'DYJetsToLL_M-50_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'DYJetsToLL_M-50_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'DYJetsToLL_M-50_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	'DYJetsToLL_M-50_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	
	'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
	'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt', 
	'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
	'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',  
	'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
	'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
	'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',  

	'WW_TuneCUETP8M1_13TeV-pythia8.txt',
	'WZ_TuneCUETP8M1_13TeV-pythia8.txt',
	'ZZ_TuneCUETP8M1_13TeV-pythia8.txt',

	'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',  
	'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',					
	'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',					
	
	# SHOULD BE AT CERN, USE eoscms.cern.ch in condor_submit.py
	# MOVE THINGS HERE TO REMIND YOURSELF WHAT IS AT CERN
	]

print '====== LJMET SUBMISSION ======'
	
relBase = os.environ['CMSSW_BASE']
print 'Relbase:',relBase

thisDir = relBase+'/src/LJMet/Com/condor_TTinclusive/' 
tarfile = relBase+'.tar'
print 'Making tar:'
if os.path.exists(tarfile):
	print 'tar already exists! Will not re-tar!'
else: 
	os.chdir(relBase)
	os.chdir('../')
	# YOU NEED TO EXCLUDE ANYTHING ELSE THAT MIGHT LIVE IN THE SAME CMSSW RELEASE
	print 'tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/LJMetSlimmer" --exclude="src/Analysis" --exclude="src/GenXsec" --exclude="src/theta" --exclude="src/tptp_2016" --exclude="tmp" -zcf'+tarfile+' '+relBase.split('/')[-1]+'/'
	os.system('tar --exclude="src/LJMet/Com/.git" --exclude="src/.git" --exclude="src/LJMetSlimmer" --exclude="src/Analysis" --exclude="src/GenXsec" --exclude="src/theta" --exclude="src/tptp_2016" --exclude="tmp" -zcf '+tarfile+' '+relBase.split('/')[-1])
	os.chdir(thisDir)

for sample in sampleList:
	os.system('python condor_submit.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsSpring16/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/2016/LJMet80X_1lep_011717 --shift '+shift)

## shift should be (one at a time): nominal, JECup, JECdown, JERup, JERdown
## If you want to use different directory names, edit lines 144 - 147 in condor_submit.py so the config is edited correctly
