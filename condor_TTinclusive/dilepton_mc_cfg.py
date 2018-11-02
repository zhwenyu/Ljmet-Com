import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import os

process = cms.Process("LJMetCom")

#Common calculator options
process.load('LJMet.Com.commonCalc_cfi')
relBase    = os.environ['CMSSW_BASE']
#condorIsMC = bool(CONDOR_ISMC)
condorIsMC = True
#relBase = str('CONDOR_RELEASE')
condorJSON = str('CONDOR_JSON')

#Dilepton calculator options
process.load('LJMet.Com.DileptonCalc_cfi')
#process.subProcesses = cms.PSet()

############################################################
#
# FWLite application options
#
process.load('LJMet.Com.ljmet_cfi')
process.ljmet = cms.PSet(
    #isMc = cms.bool(True),
    isMc = cms.bool(condorIsMC),
    keepFullMChistory = cms.bool(condorIsMC),
    runs = cms.vint32([])
)

process.ljmet.excluded_calculators = cms.vstring(
	'JetSubCalc',
	'BestCalc', #
	'DeepAK8Calc', #
	'PileUpCalc',
	'LjetsTopoCalcMinPz', #
        'LjetsTopoCalcNew' #
	'StopCalc',
	'PdfCalc',
	'ChargedHiggsCalc',
	'TprimeCalc',
	'TpTpCalc', #
	'LjetsTopoCalc',
	'WprimeCalc',
        'WprimeBoostedCalc',
	'singleLepCalc',
        'BTagSFCalc',
        'CATopoCalc',
	#'bjornCalc',
	) 

# Jet substructure calculator options -- see src/JetSubCalc.cc
#process.load('LJMet.Com.JetSubCalc_cfi')
#process.JetSubCalc.useHTT = cms.bool(False)
#process.JetSubCalc.killHF = cms.bool(False)
#process.JetSubCalc.doNewJEC = cms.bool(True)
#process.JetSubCalc.useL2L3Mass = cms.bool(True)
#process.JetSubCalc.isMc = process.ljmet.isMc
#process.JetSubCalc.JECup = cms.bool(False)
#process.JetSubCalc.JECdown = cms.bool(False)
#process.JetSubCalc.MCL2JetParAK8 = cms.string('../data/Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt')
#process.JetSubCalc.MCL3JetParAK8 = cms.string('../data/Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt')
#process.JetSubCalc.DataL2JetParAK8 = cms.string('../data/Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt')
#process.JetSubCalc.DataL3JetParAK8 = cms.string('../data/Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt')
#process.JetSubCalc.DataL2L3JetParAK8 = cms.string('../data/Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt')
#process.JetSubCalc.UncertaintyAK8 = cms.string('../data/Spring16_25nsV6_DATA_Uncertainty_AK8PFchs.txt')

#process.DileptonCalc.isMc     = process.ljmet.isMc
process.DileptonCalc.isMc = cms.bool(condorIsMC)
process.DileptonCalc.keepFullMChistory = cms.bool(condorIsMC)
process.DileptonCalc.dataType = cms.string('ElMu')
process.DileptonCalc.UseElMVA = cms.bool(True)

############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('DileptonSelector'),
    #selection = cms.string('singleLepSelector'),
    isMc              = cms.bool(True),
    keepFullMChistory = cms.bool(True),
    debug             = cms.bool(False),
    #debug             = cms.bool(True),
    # cuts
    #HLT
    trigger_cut              = cms.bool(False),
    dump_trigger             = cms.bool(True),

    #Can use same trigger paths for data and MC since MC is always one of the data versions
    trigger_path_ee          = cms.vstring('HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1', 'HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v1'),
    #trigger_path_ee          = cms.vstring('HLT_Ele27_Ele37_CaloIdL_MW_v1')
    
    trigger_path_em          = cms.vstring('HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_IsoVL_v1', 'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_IsoVL_v1'),

    trigger_path_mm          = cms.vstring('HLT_Mu17_Mu8_v1', 'HLT_Mu17_TkMu8_v1', 'HLT_Mu30_TkMu11_v1', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1'),    


    ElMVAweightFiles = cms.vstring(
        relBase+'/src/LJMet/Com/weights/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EE_10_2017_puinfo_BDT.weights.xml',
        ),

    ElMVAweightFiles_iso = cms.vstring(
        relBase+'/src/LJMet/Com/weights/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml',
        ),


    pv_cut                   = cms.bool(True),
    metfilters		     = cms.bool(True), #
    hbhe_cut                 = cms.bool(True),
    hbhe_cut_value           = cms.string("Run2Loose"),
    hbheiso_cut              = cms.bool(True),
    eesc_cut                 = cms.bool(True),
    ecalTP_cut               = cms.bool(True),
    goodVtx_cut              = cms.bool(True),
    badMuon_cut              = cms.bool(True),
    badChargedHadron_cut     = cms.bool(True),    			     
    flag_tag                 = cms.InputTag('TriggerResults::PAT'),
    cscHalo_cut              = cms.bool(True),
    jet_cuts                 = cms.bool(False),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(5.0),
    min_jet                  = cms.int32(4),
    max_jet                  = cms.int32(4000),

    muon_cuts                = cms.bool(True),
    min_muon                 = cms.int32(0),
    muon_minpt               = cms.double(10.0),
    muon_maxeta              = cms.double(2.4),
    max_muon                 = cms.int32(10),

    electron_cuts            = cms.bool(True),
    min_electron             = cms.int32(0),
    electron_minpt           = cms.double(10.0),
    electron_maxeta          = cms.double(2.4),
    max_electron             = cms.int32(10),
    UseElMVA                 = cms.bool(True),

    min_lepton               = cms.int32(2),

    #B-tagging cuts
    btagOP                   = cms.string('MEDIUM'),
    bdisc_min                = cms.double(0.4941),
    DeepCSVfile              = cms.string(relBase+'/src/LJMet/Com/data/DeepCSV_94XSF_V3_B_F.csv'),
    DeepCSVSubjetfile        = cms.string(relBase+'/src/LJMet/Com/data/subjet_DeepCSV_94XSF_V3_B_F.csv'),
    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(0.0),
    btag_cuts                = cms.bool(False),
    #btagOP                  = cms.string("CSVM"),
    btag_1                   = cms.bool(True),
    btag_2                   = cms.bool(True),
    btag_3                   = cms.bool(False),

    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    jet_collection           = cms.InputTag('slimmedJets'),
    muon_collection          = cms.InputTag('slimmedMuons'),
    electron_collection      = cms.InputTag('slimmedElectrons'),
    met_collection           = cms.InputTag('slimmedMETs'),

    #This is for uncertainty calculations   
    JEC_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt'),
    JERSF_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),
    JER_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),
    JERAK8_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt'),


    JECup		     = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),

    #new jet energy corrections
    doNewJEC                 = cms.bool(True),
    doLepJetCleaning         = cms.bool(True),

    MCL1JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017_V6_MC_L1FastJet_AK4PFchs.txt'),
    MCL2JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017_V6_MC_L2Relative_AK4PFchs.txt'),
    MCL3JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017_V6_MC_L3Absolute_AK4PFchs.txt'),

    MCL1JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017_V6_MC_L1FastJet_AK8PFPuppi.txt'),
    MCL2JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017_V6_MC_L2Relative_AK8PFPuppi.txt'),
    MCL3JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017_V6_MC_L3Absolute_AK8PFPuppi.txt'),

    DataL1JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK4PFchs.txt'),
    DataL2JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L2Relative_AK4PFchs.txt'),
    DataL3JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK4PFchs.txt'),
    DataResJetPar            = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK4PFchs.txt'),

    DataL1JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK8PFPuppi.txt'),
    DataL2JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L2Relative_AK8PFPuppi.txt'),
    DataL3JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK8PFPuppi.txt'),
    DataResJetParAK8         = cms.string(relBase+'/src/LJMet/Com/data/Fall17V6/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK8PFPuppi.txt'),


    #triggerstudy info
    doTriggerStudy           = cms.bool(True),
    #doTriggerStudy           = cms.bool(False),
    TriggerBits              = cms.InputTag("TriggerResults","","HLT"),
    TriggerObjects           = cms.InputTag("selectedPatTrigger"),


)

#######################################################
#
# Input files
#


process.inputs = cms.PSet (
   nEvents    = cms.int32(10000),
   #nEvents    = cms.int32(-1),
   skipEvents = cms.int32(0),	
   lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
       fileNames  = cms.vstring(
	#INFILES #i think this was an old variable name
	#CONDOR_FILELIST
	#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/TTTT_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/026CBACC-4553-E811-A711-0CC47A13CCEE.root', 
	'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/X53X53To2L2Nu_M-1200_LH_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/2AF651B9-5A44-E811-8202-0CC47A7C345C.root',
     ) 
)



# JSON
if (not process.ljmet.isMc==cms.bool(True)):
    JsonFile = JSONFILE
    JsonFile = relBase+'/src/LJMet/Com/data/json/'+condorJSON
    myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
    process.inputs.lumisToProcess.extend(myList)
       
        
        
#######################################################
#
# Output
#

import os
try: scratchDir=os.environ['/uscms_data/d3/bburkle/working_area/CMSSW_9_4_6_patch1/src/LJMet/Com/mc_out/']
except: scratchDir='.'
if not os.path.isdir(scratchDir): os.system("mkdir -p "+scratchDir)
process.outputs = cms.PSet (
    #outputName = cms.string(scratchDir+'/PREFIX_JOBID'),
    outputName = cms.string('ljmet_X53X53LH-m1200_2017'),
    #outputName = cms.string('CONDOR_OUTFILE'),
    treeName   = cms.string('ljmet'),
)

#######################################################
#
# Object selector options
#

# Primary vertex
process.load('PhysicsTools.SelectorUtils.pvSelector_cfi')
process.pvSelector.pvSrc   = cms.InputTag('offlineSlimmedPrimaryVertices')
process.pvSelector.minNdof = cms.double(4.0)
process.pvSelector.maxZ    = cms.double(24.0)
process.pvSelector.maxRho  = cms.double(2.0)

# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')
