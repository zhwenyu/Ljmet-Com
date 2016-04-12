import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

# Define the base process
process = cms.Process("LJMetCom")

#Arguments from condor submit script which are used more than once
condorIsMC = bool(True)
relBase    = str('/uscms_data/d3/jmanagan/CleanLJMet76X')
condorJSON = str('CONDOR_JSON')
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

# common calculator options -- see src/CommonCalc.cc
process.load('LJMet.Com.commonCalc_cfi')

# singleLep calculator options -- see src/singleLepCalc.cc
process.load('LJMet.Com.singleLepCalc_cfi')
process.singleLepCalc.isMc              = cms.bool(condorIsMC)
process.singleLepCalc.keepFullMChistory = cms.bool(condorIsMC)
process.singleLepCalc.UseElMVA          = cms.bool(True)
process.singleLepCalc.saveLooseLeps     = cms.bool(False)

# Jet substructure calculator options -- see src/JetSubCalc.cc
process.load('LJMet.Com.JetSubCalc_cfi')
process.JetSubCalc.useHTT = cms.bool(False)
process.JetSubCalc.killHF = cms.bool(False)
process.JetSubCalc.doNewJEC = cms.bool(True)
process.JetSubCalc.useL2L3Mass = cms.bool(True)
process.JetSubCalc.isMc = cms.bool(condorIsMC)
process.JetSubCalc.JECup = cms.bool(False)
process.JetSubCalc.JECdown = cms.bool(False)
process.JetSubCalc.MCL2JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L2Relative_AK8PFchs.txt')
process.JetSubCalc.MCL3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L3Absolute_AK8PFchs.txt')
process.JetSubCalc.DataL2JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L2Relative_AK8PFchs.txt')
process.JetSubCalc.DataL3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L3Absolute_AK8PFchs.txt')
process.JetSubCalc.DataL2L3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt')
process.JetSubCalc.UncertaintyAK8 = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_Uncertainty_AK8PFchs.txt')

############################################################
#
# Event selector options -- see src/singleLepEventSelector.cc
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
        'HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v3',
        'HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v3',
        'HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3',
        'HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v2',
        'HLT_Ele15_IsoVVVL_PFHT600_v3',
        'HLT_Ele15_IsoVVVL_PFHT350_v2',
        'HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2',
        'HLT_Ele22_eta2p1_WPLoose_Gsf_v3',
        'HLT_Ele22_eta2p1_WPTight_Gsf_v3',
        'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v3',
        'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3',
        'HLT_Ele23_WPLoose_Gsf_v3',
        'HLT_Ele23_WPLoose_Gsf_TriCentralPFJet50_40_30_v2',
        'HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v3',
        'HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v2',
        'HLT_Ele27_WPLoose_Gsf_v1',
        'HLT_Ele27_WPLoose_Gsf_TriCentralPFJet50_40_30_v1',
        'HLT_Ele27_eta2p1_WPLoose_Gsf_v2',
        'HLT_Ele27_eta2p1_WPTight_Gsf_v2',
        'HLT_Ele32_eta2p1_WPTight_Gsf_v2',
        'HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v3',
        'HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v3',
        'HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_v1',
        'HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v3',
        'HLT_Ele105_CaloIdVT_GsfTrkIdT_v3',
        'HLT_Ele115_CaloIdVT_GsfTrkIdT_v2',
        ),
    mctrigger_path_mu = cms.vstring(
        'HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v2',
        'HLT_Mu15_IsoVVVL_PFHT600_v3',
        'HLT_Mu15_IsoVVVL_PFHT350_v2',
        'HLT_Mu16_eta2p1_MET30_v1',
        'HLT_IsoMu16_eta2p1_MET30_v1',
        'HLT_Mu17_TrkIsoVVL_v2',
        'HLT_Mu17_v2',
        'HLT_IsoMu17_eta2p1_v3',
        'HLT_IsoMu18_v2',
        'HLT_IsoTkMu18_v2',
        'HLT_OldIsoTkMu18_v2',
        'HLT_OldIsoMu18_v1',
        'HLT_IsoMu18_TriCentralPFJet50_40_30_v2',
        'HLT_Mu20_v2',
        'HLT_TkMu20_v2',
        'HLT_IsoMu20_v3',
        'HLT_IsoTkMu20_v4',
        'HLT_IsoMu22_v2',
        'HLT_IsoMu22_TriCentralPFJet50_40_30_v2',
        'HLT_IsoTkMu22_v2',
        'HLT_Mu24_eta2p1_v2',
        'HLT_TkMu24_eta2p1_v2',
        'HLT_IsoMu27_v3',
        'HLT_Mu27_v2',
        'HLT_TkMu27_v2',
        'HLT_IsoTkMu27_v3',
        'HLT_Mu27_TkMu8_v2',
        'HLT_Mu30_TkMu11_v2',
        'HLT_Mu30_eta2p1_PFJet150_PFJet50_v1',
        'HLT_Mu40_TkMu11_v2',
        'HLT_Mu40_eta2p1_PFJet200_PFJet50_v3',
        'HLT_Mu45_eta2p1_v2',
        'HLT_Mu50_v2',
        'HLT_Mu55_v1',
        'HLT_Mu300_v1',
        'HLT_Mu350_v1',
        ),
    
    trigger_path_el = cms.vstring(''),
    trigger_path_mu = cms.vstring(''),   
    
    # PV cuts
    pv_cut         = cms.bool(True),
    hbhe_cut       = cms.bool(True),
    hbheiso_cut    = cms.bool(True),
    hbhe_cut_value = cms.string('Run2Loose'),
    csc_cut        = cms.bool(False),
    flag_tag       = cms.InputTag('TriggerResults::PAT'),
    eesc_cut       = cms.bool(True),
    
    # Jet cuts
    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(5.0),
    min_jet                  = cms.int32(2),
    max_jet                  = cms.int32(4000),
    leading_jet_pt           = cms.double(30.0),

    # muon cuts
    muon_cuts                = cms.bool(True),
    min_muon                 = cms.int32(0),
    # use a code based selector? Right now "false" for all
    muon_selector            = cms.bool(False),
    muon_selector_medium     = cms.bool(False),
    loose_muon_selector      = cms.bool(False),
    loose_muon_selector_tight = cms.bool(False),
    # turn on miniIso and set cuts
    muon_useMiniIso          = cms.bool(True),
    muon_miniIso             = cms.double(0.2),
    loose_muon_miniIso       = cms.double(0.4),
    # or use relative isolation
    muon_reliso              = cms.double(0.2),
    loose_muon_reliso        = cms.double(0.4),
    # choose min/max pt and eta
    muon_minpt               = cms.double(30.0),
    muon_maxeta              = cms.double(2.4),
    loose_muon_minpt         = cms.double(10.0),
    loose_muon_maxeta        = cms.double(2.4),

    # electron cuts
    electron_cuts            = cms.bool(True),
    min_electron             = cms.int32(0),
    # turn on MVA identification, or default to cut-based
    UseElMVA                 = cms.bool(True),
    tight_electron_mva_cuts  = cms.vdouble(0.967083,0.929117,0.726311), # ~80% el efficiency WP
    loose_electron_mva_cuts  = cms.vdouble(0.913286,0.805013,0.358969), # ~90% el efficiency WP
    # with MVA ID we use mini-isolation, set the cut
    electron_miniIso         = cms.double(0.1),
    loose_electron_miniIso   = cms.double(0.4),
    # choose min/max pt and eta
    electron_minpt           = cms.double(30.0),
    electron_maxeta          = cms.double(2.1),
    loose_electron_minpt     = cms.double(10.0),
    loose_electron_maxeta    = cms.double(2.1),

    ElMVAweightFiles = cms.vstring(
        relBase+'/src/LJMet/Com/weights/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml',
        relBase+'/src/LJMet/Com/weights/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml',
        ),

    # more lepton cuts
    min_lepton               = cms.int32(1),    # checks (N tight mu + N tight el) >= cut
    max_lepton               = cms.int32(1),    # checks (N tight mu + N tight el) <= cut
    min_loose_lepton         = cms.int32(1),    # checks (N loose mu + N loose el) >= cut
    max_loose_lepton         = cms.int32(1),    # checks (N loose mu + N loose el) <= cur
    second_lepton_veto       = cms.bool(True),  # checks (N tight lep > 0) AND (N loose lep > 0), vetoes if there are loose leptons.
    tau_veto		     = cms.bool(False),
    
    # MET cuts
    met_cuts                 = cms.bool(True),
    min_met                  = cms.double(20.0),
    max_met                  = cms.double(100000.0),
    
    # Btagging cuts
    btagOP                   = cms.string('CSVM'),
    btag_min_discr           = cms.double(0.800),
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
    
    # Jet corrections are read from txt files
    JECup                    = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    doNewJEC                 = cms.bool(True),
    doLepJetCleaning         = cms.bool(True),
    LepJetDR                 = cms.double(0.4),
    CleanLooseLeptons        = cms.bool(False),

    # As of Fall15_25nsV2, JEC uncertainty and JERSF are identical for AK4 & AK8 PFchs
    JEC_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt'),
    JERSF_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
    JER_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
    JERAK8_txtfile = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_MC_PtResolution_AK8PFchs.txt'),
    
    MCL1JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L1FastJet_AK4PFchs.txt'),
    MCL2JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L2Relative_AK4PFchs.txt'),
    MCL3JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L3Absolute_AK4PFchs.txt'),

    MCL1JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L1FastJet_AK8PFchs.txt'),
    MCL2JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L2Relative_AK8PFchs.txt'),
    MCL3JetParAK8            = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV1_MC_L3Absolute_AK8PFchs.txt'),

    DataL1JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt'),
    DataL2JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt'),
    DataL3JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt'),
    DataResJetPar            = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt'),

    DataL1JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L1FastJet_AK8PFchs.txt'),
    DataL2JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L2Relative_AK8PFchs.txt'),
    DataL3JetParAK8          = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L3Absolute_AK8PFchs.txt'),
    DataResJetParAK8         = cms.string(relBase+'/src/LJMet/Com/data/Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt')

    )


#######################################################
#
# Input files
#

process.inputs = cms.PSet (
    nEvents    = cms.int32(5000),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    fileNames  = cms.vstring(
        #'root://cmseos.fnal.gov//store/mc/RunIIFall15MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/FEACEF28-8EBC-E511-A9E2-001E67398683.root',
        #'root://cmseos.fnal.gov//store/mc/RunIIFall15MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/FA364BB4-17BC-E511-A3D9-008CFA197DF8.root'
        #'root://cmseos.fnal.gov//store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/FC27F1DE-4DC2-E511-8675-000F530E4774.root',
        'root://cmseos.fnal.gov//store/mc/RunIIFall15MiniAODv2/TprimeTprime_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/F06A5543-0CB9-E511-A3B4-B083FECF8ACE.root',
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
    outputName = cms.string('testmc'),
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

# Tight muon -- not used
process.load('LJMet.Com.pfMuonSelector_cfi') 

# Loose muon -- not used
process.LoosepfMuonSelector = process.pfMuonSelector.clone()

# Tight electron for 25ns --- only used if MVA ID is turned off
process.load('LJMet.Com.TopElectronSelector_cfi')
process.TopElectronSelector.version = cms.string('NONE')
process.TopElectronSelector.deta_EB	  = cms.double(0.00926)
process.TopElectronSelector.dphi_EB	  = cms.double(0.0336)
process.TopElectronSelector.sihih_EB	  = cms.double(0.0101)
process.TopElectronSelector.hoe_EB	  = cms.double(0.0597)
process.TopElectronSelector.d0_EB	  = cms.double(0.0111)
process.TopElectronSelector.dZ_EB	  = cms.double(0.0466)
process.TopElectronSelector.ooemoop_EB	  = cms.double(0.012)
process.TopElectronSelector.reliso_EB	  = cms.double(0.0354)
process.TopElectronSelector.deta_EE	  = cms.double(0.00724)
process.TopElectronSelector.dphi_EE	  = cms.double(0.0918)
process.TopElectronSelector.sihih_EE	  = cms.double(0.0279)
process.TopElectronSelector.hoe_EE	  = cms.double(0.0615)
process.TopElectronSelector.d0_EE	  = cms.double(0.0351)
process.TopElectronSelector.dZ_EE	  = cms.double(0.417)
process.TopElectronSelector.ooemoop_EE	  = cms.double(0.00999)
process.TopElectronSelector.reliso_EE	  = cms.double(0.0646)
                   	       
#Loose electron for 25ns --- only used if MVA ID is turned off 
process.LooseTopElectronSelector = process.TopElectronSelector.clone()
process.LooseTopElectronSelector.version = cms.string('NONE')
process.LooseTopElectronSelector.deta_EB	  = cms.double(0.0105)
process.LooseTopElectronSelector.dphi_EB	  = cms.double(0.115)
process.LooseTopElectronSelector.sihih_EB	  = cms.double(0.0103)
process.LooseTopElectronSelector.hoe_EB           = cms.double(0.104)
process.LooseTopElectronSelector.d0_EB		  = cms.double(0.0261)
process.LooseTopElectronSelector.dZ_EB		  = cms.double(0.41) 
process.LooseTopElectronSelector.ooemoop_EB	  = cms.double(0.102)
process.LooseTopElectronSelector.reliso_EB	  = cms.double(0.0893)
process.LooseTopElectronSelector.deta_EE	  = cms.double(0.00814)
process.LooseTopElectronSelector.dphi_EE	  = cms.double(0.182)
process.LooseTopElectronSelector.sihih_EE	  = cms.double(0.0301)
process.LooseTopElectronSelector.hoe_EE	          = cms.double(0.0897)
process.LooseTopElectronSelector.d0_EE		  = cms.double(0.118)
process.LooseTopElectronSelector.dZ_EE		  = cms.double(0.822)
process.LooseTopElectronSelector.ooemoop_EE	  = cms.double(0.126)
process.LooseTopElectronSelector.reliso_EE	  = cms.double(0.121)





