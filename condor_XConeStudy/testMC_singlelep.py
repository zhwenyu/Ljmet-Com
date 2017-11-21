import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import os

# Define the base process
process = cms.Process("LJMetCom")

#Arguments from condor submit script which are used more than once
condorIsMC = bool(True)
relBase    = os.environ['CMSSW_BASE']
condorJSON = str('none')
############################################################
#
# FWLite application options
process.load('LJMet.Com.ljmet_cfi')
process.ljmet.isMc = cms.bool(condorIsMC)

# Exclude some unnecessary calculators from the process
process.ljmet.excluded_calculators = cms.vstring(
    'PileUpCalc',
    'BTagSFCalc',
    'TprimeCalc',
    'CATopoCalc',
    'DileptonCalc',
    'DileptonEventSelector',
    'StopCalc',
    'PdfCalc',
    'ChargedHiggsCalc',
    'LjetsTopoCalc',
    'WprimeCalc'
    ) 

# common calculator options
process.load('LJMet.Com.commonCalc_cfi')

# singleLep calculator options
process.load('LJMet.Com.singleLepCalc_cfi')
process.singleLepCalc.isMc              = cms.bool(condorIsMC)
process.singleLepCalc.keepFullMChistory = cms.bool(condorIsMC)
process.singleLepCalc.UseElMVA          = cms.bool(True)
process.singleLepCalc.saveLooseLeps     = cms.bool(True)
process.singleLepCalc.saveGenHT         = cms.bool(False)

# Jet substructure calculator options
process.load('LJMet.Com.JetSubCalc_cfi')
process.JetSubCalc.killHF = cms.bool(False)
process.JetSubCalc.doNewJEC = cms.bool(True)
process.JetSubCalc.useL2L3Mass = cms.bool(True)
process.JetSubCalc.isMc = cms.bool(condorIsMC)
process.JetSubCalc.MCL2JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt')
process.JetSubCalc.MCL3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt')
process.JetSubCalc.MCSF = cms.string(relBase+'/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_SF_AK4PFchs.txt')
process.JetSubCalc.DataL2JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFchs.txt')
process.JetSubCalc.DataL3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFchs.txt')
process.JetSubCalc.DataL2L3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFchs.txt')
process.JetSubCalc.UncertaintyAK8 = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_Uncertainty_AK8PFchs.txt')


# XCone calculator options
process.load('LJMet.Com.XConeCalc_cfi')
process.XConeCalc.XConeR			 = cms.double(0.4)
# process.XConeCalc.XConeNumJets		 = cms.int32(10) ##ttbar full hadronic : 6 , TT (with fully hadronic top): 10 ! But it doesn't matter if youa re using VarNumJets
process.XConeCalc.XConeBeta			 = cms.double(2.0)
process.XConeCalc.usePFchs			 = cms.bool(True)
process.XConeCalc.doPUPPI			 = cms.bool(False)
process.XConeCalc.isMc				 = cms.bool(condorIsMC)
process.XConeCalc.doGenXCone		 = cms.bool(False)
process.XConeCalc.DEBUG			     = cms.bool(True)
process.XConeCalc.saveJetConst		 = cms.bool(False)
process.XConeCalc.MCL1JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt')
process.XConeCalc.MCL2JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt')
process.XConeCalc.MCL3JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt')
process.XConeCalc.DataL1JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt')
process.XConeCalc.DataL2JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt')
process.XConeCalc.DataL3JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt')
process.XConeCalc.DataResJetPar            = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt')

############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(
    
    selection = cms.string('singleLepSelector'),
    
    # Define cuts -- variable names are strings searched by src/singleLepEventSelector.cc
    
    debug  = cms.bool(False),
    
    isMc  = cms.bool(condorIsMC),
    keepFullMChistory = cms.bool(condorIsMC),
    doLaserCalFilt  = cms.bool(False),
    
    # Trigger cuts
    trigger_cut  = cms.bool(True),
    dump_trigger = cms.bool(False),

    mctrigger_path_el = cms.vstring(
		#from Julies singleLep
        'digitisation_step',
        'HLT_Ele15_IsoVVVL_PFHT350_v',            
        'HLT_Ele15_IsoVVVL_PFHT400_v',            
        'HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v',   
        'HLT_Ele32_eta2p1_WPTight_Gsf_v',         
        'HLT_Ele30_WPTight_Gsf_v',                
        'HLT_Ele32_WPTight_Gsf_v',                
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_v',
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v',
		#from Clints fakerate study
		'HLT_Ele17_CaloIdL_GsfTrkIdVL_v',
        ),
    mctrigger_path_mu = cms.vstring(
		#from Julies singleLep
        'digitisation_step',
        'HLT_IsoMu24_v',                     
        'HLT_IsoTkMu24_v',                   
        'HLT_Mu50_v',                        
        'HLT_TkMu50_v',                      
        'HLT_Mu15_IsoVVVL_PFHT350_v',        
        'HLT_Mu15_IsoVVVL_PFHT400_v',        
        'HLT_Mu40_eta2p1_PFJet200_PFJet50_v',
		#from Clints fakerate study
        'HLT_Mu17_v',
        ),
    
    trigger_path_el = cms.vstring(''),
    trigger_path_mu = cms.vstring(''),   
    
    # PV cuts
    pv_cut         = cms.bool(True),
    flag_tag       = cms.InputTag('TriggerResults::PAT'),
    metfilters     = cms.bool(True),
    
    # Jet cuts
    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(2.5),
    min_jet                  = cms.int32(1),
    max_jet                  = cms.int32(4000),
    leading_jet_pt           = cms.double(30.0),

    # muon cuts
    muon_cuts                = cms.bool(True),
    min_muon                 = cms.int32(0),
    muon_minpt               = cms.double(20.0),
    muon_maxeta              = cms.double(2.4),
    muon_useMiniIso          = cms.bool(True),
    muon_miniIso             = cms.double(0.2),
    loose_muon_miniIso       = cms.double(0.4),
    loose_muon_minpt         = cms.double(20.0),
    loose_muon_maxeta        = cms.double(2.4),
    # specify IP cuts
    muon_dxy                 = cms.double(0.2),
    muon_dz                  = cms.double(0.5),
    loose_muon_dxy           = cms.double(999999.),
    loose_muon_dz            = cms.double(999999.),
 
    # electron cuts
    electron_cuts            = cms.bool(True),
    min_electron             = cms.int32(0),
    electron_minpt           = cms.double(20.0),
    electron_maxeta          = cms.double(2.4),
    electron_useMiniIso      = cms.bool(True),
    electron_miniIso         = cms.double(0.1),
    loose_electron_miniIso   = cms.double(0.4),
    loose_electron_minpt     = cms.double(20.0),
    loose_electron_maxeta    = cms.double(2.4),
    UseElMVA                 = cms.bool(True),
    tight_electron_mva_cuts  = cms.vdouble(0.674,0.744,0.170), # 80X WP80 to recover efficiency of 74X WP80
    loose_electron_mva_cuts  = cms.vdouble(-0.041,0.383,-0.515), # 80X WP90 to recover efficiency of 74X WP90
    #tight_electron_mva_cuts  = cms.vdouble(0.967083,0.929117,0.726311), # ~80% el efficiency WP 74X
    #loose_electron_mva_cuts  = cms.vdouble(0.913286,0.805013,0.358969), # ~90% el efficiency WP 74X

    ElMVAweightFiles = cms.vstring(
        relBase+'/src/LJMet/Com/weights/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml',
        ),

    ElMVAweightFiles_alt = cms.vstring(
        relBase+'/src/LJMet/Com/weights/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml',
        relBase+'/src/LJMet/Com/weights/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml',
        relBase+'/src/LJMet/Com/weights/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml',
        ),

    # more lepton cuts
    min_lepton               = cms.int32(1),    # checks (N tight mu + N tight el) >= cut
    max_lepton               = cms.int32(1),    # checks (N tight mu + N tight el) <= cut
    min_loose_lepton         = cms.int32(1),
    max_loose_lepton         = cms.int32(1),
    second_lepton_veto       = cms.bool(False),  #checks (N tight lep > 0) AND (N loose lep > 0), vetoes if there are loose leptons.
    tau_veto		     = cms.bool(True),
    
    # MET cuts
    met_cuts                 = cms.bool(True),
    min_met                  = cms.double(20.0),
    max_met                  = cms.double(99999999999.0),
    
    # Btagging cuts
    btagOP                   = cms.string('CSVM'),
    btag_min_discr           = cms.double(0.8484),
    btag_cuts                = cms.bool(False),
    btag_1                   = cms.bool(False),
    btag_2                   = cms.bool(False),
    btag_3                   = cms.bool(False),
    
    # Define the branch names of object collections in the input miniAOD file
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    jet_collection           = cms.InputTag('slimmedJets'),
    muon_collection          = cms.InputTag('slimmedMuons'),
    electron_collection      = cms.InputTag('slimmedElectrons'),
    tau_collection	     = cms.InputTag('slimmedTaus'),
    met_collection           = cms.InputTag('slimmedMETs'),
    
    # Jet corrections are read from txt files which need updating!
    JEC_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt'),
    JERSF_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),
    JER_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),
    JERAK8_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt'),

    JECup                    = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    doNewJEC                 = cms.bool(True),
    doLepJetCleaning         = cms.bool(True),
    CleanLooseLeptons        = cms.bool(True),
    LepJetDR                 = cms.double(0.4),
    
    MCL1JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt'),
    MCL2JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt'),
    MCL3JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt'),

    MCL1JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt'),
    MCL2JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt'),
    MCL3JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt'),

    DataL1JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt'),
    DataL2JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt'),
    DataL3JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt'),
    DataResJetPar            = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt'),

    DataL1JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK8PFchs.txt'),
    DataL2JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK8PFchs.txt'),
    DataL3JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK8PFchs.txt'),
    DataResJetParAK8         = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK8PFchs.txt'),


    # Unused parameters
    muon_reliso              = cms.double(0.2),
    muon_selector            = cms.bool(False),
    muon_selector_medium     = cms.bool(False),
    loose_muon_selector      = cms.bool(False),
    loose_muon_selector_tight = cms.bool(False),
    loose_muon_reliso        = cms.double(0.4),
    electron_CutsPlusMVA     = cms.bool(False),
    BTagUncertUp             = cms.bool(False), # no longer needed
    BTagUncertDown           = cms.bool(False), # no longer needed
 
    )


#######################################################
#
# Input files
#

process.inputs = cms.PSet (
#     nEvents    = cms.int32(4000),
    nEvents    = cms.int32(1000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring(
# 		'root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/029929BB-ACBD-E611-96C1-002590DE6C9A.root',
		'root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
        )
    )


# JSON
if (not process.ljmet.isMc==cms.bool(True)):
    JsonFile = relBase+'/src/LJMet/Com/data/json/'+condorJSON
    myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
    process.inputs.lumisToProcess.extend(myList)
    
    
#######################################################
#
# Output
#
process.outputs = cms.PSet (
#     outputName = cms.string('testmc_LJMet'),
#     outputName = cms.string('QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_testMC_singlelep_LJMet'),
    outputName = cms.string('TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_testMC_singlelep_LJMetXCone'),
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

# Tight muon
process.load('LJMet.Com.pfMuonSelector_cfi') #not used

# Loose muon
process.LoosepfMuonSelector = process.pfMuonSelector.clone() #not used

# Tight electron for 25ns
process.load('LJMet.Com.TopElectronSelector_cfi')
process.TopElectronSelector.version = cms.string('NONE') #not used
                   	       
#Loose electron for 25ns
process.LooseTopElectronSelector = process.TopElectronSelector.clone()
process.LooseTopElectronSelector.version = cms.string('NONE') #not used
