import os,sys,datetime

shift = sys.argv[1]

sampleList=[
	'TprimeTprime_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'TprimeTprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',

	'BprimeBprime_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1100_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',
	'BprimeBprime_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',

	'X53X53_M-700_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',                               
	'X53X53_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',				      
	'X53X53_M-900_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',				      
	'X53X53_M-1000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1100_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1300_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1500_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-700_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',				      
	'X53X53_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',				      
	'X53X53_M-900_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',				      
	'X53X53_M-1000_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1100_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1200_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1300_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1500_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
	'X53X53_M-1600_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_25ns.txt',			      
												      
#	'TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns.txt',				      
#	'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',				      
	'TT_TuneCUETP8M1_13TeV-powheg-pythia8_25ns.txt',					      
	'TT_TuneCUETP8M1_13TeV-powheg-pythia8_highstats_25ns.txt',				      
	'TT_Mtt-1000toInf_TuneCUETP8M1_13TeV-powheg-pythia8_25ns.txt',				      
	'TT_Mtt-700to1000_TuneCUETP8M1_13TeV-powheg-pythia8_25ns.txt',				      
#	'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns.txt',				      
#	'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',				      
	'WW_TuneCUETP8M1_13TeV-pythia8_25ns.txt',						      
	'WZ_TuneCUETP8M1_13TeV-pythia8_25ns.txt',						      
	'ZZ_TuneCUETP8M1_13TeV-pythia8_25ns.txt',						      
#	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_25ns.txt',			      
	'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_25ns.txt',		      
	'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_25ns.txt',		      
	'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns.txt',		      
	'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns.txt',		      
	'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_25ns.txt',		      
	'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_25ns.txt',			      
	'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_25ns.txt',					      
	'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_25ns.txt',			      
												      
	'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',			      
	'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',              	      
												      
	'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',		      
	'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',		      
	'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',		      
	'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',		      
	'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',		      
	'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',		      
	'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns.txt',		      
												      
#	'TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8_25ns.txt',				      
#	'TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8_25ns.txt',				      
#	'ST_tW_top_5f_scaleup_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns.txt',	      
#	'ST_tW_top_5f_scaledown_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns.txt',	      
#	'ST_tW_antitop_5f_scaleup_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns.txt',	      
#	'ST_tW_antitop_5f_scaledown_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_25ns.txt',      

]

relBase = os.environ['CMSSW_BASE']
thisDir = relBase+'/src/LJMet/Com/condor/' 
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
	os.system('python TTsinglelep_condorsubmit.py --useMC True --sample '+sample.split('.')[0]+' --fileList '+thisDir+'fileLists/'+sample+' --submit True --inputTar '+tarfile+' --outDir /eos/uscms/store/user/lpcljm/LJMet_1lepTT_020516 --shift '+shift)
