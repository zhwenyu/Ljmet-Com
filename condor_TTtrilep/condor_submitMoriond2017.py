import os,sys,datetime

shift = sys.argv[1]

sampleList=[

	###SIGNAL:
	'TprimeTprime_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  
	'TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  
	'TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8.txt',  								       
	'TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 
	'TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8.txt', 

	###BKG:
	'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.txt',				  
										  
# 	'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.txt',	  
# 	'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',	  
# 	'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',	  
# 	'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt', 
# 	'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.txt',	  
	
	'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',
# 	'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
# 	
	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',
# 	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	
# 	'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt',   
# 	'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
# 	'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
# 	'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt', 
# 	'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
# 	'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
# 	'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',   
# 	'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',  

	'WW_TuneCUETP8M1_13TeV-pythia8.txt',
	'WZ_TuneCUETP8M1_13TeV-pythia8.txt',
	'ZZ_TuneCUETP8M1_13TeV-pythia8.txt',

	#'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',  
	#'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	#'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	#'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	#'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',	
	#'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt',
	
	#Additionally for multilepton:
	'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.txt',
	'ZZTo4L_13TeV_powheg_pythia8.txt',
	'WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',
	'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',       
						
	'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.txt',       
	'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.txt',	  
	'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',			  
	'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt',			  

	]

# if shift == 'nominal':
# 	# USE cmsxrootd.fnal.gov for scaledown (or both)
# 	sampleList.append('TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8.txt'),
# 	sampleList.append('TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8.txt'),

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
outdir = 'LJMet80x_3lepTT_Full2016_Moriond17_newJEC_newElMVA_'+date+'_rizki'

for sample in sampleList:
	os.system('python condor_submit.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileListsMoriond17/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/'+outdir+' --shift '+shift)

## shift should be (one at a time): nominal, JECup, JECdown, JERup, JERdown
## If you want to use different directory names, edit lines 144 - 147 in condor_submit.py so the config is edited correctly
