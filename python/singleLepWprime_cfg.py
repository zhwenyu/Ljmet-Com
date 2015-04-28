import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("LJMetCom")

relBase    = str('/data1/speer/tblsm/cmssw/git-transition/new-git_5_3_6/')
############################################################
#
# FWLite application options
process.load('LJMet.Com.ljmet_cfi')
process.ljmet.isMc = cms.bool(True)
process.ljmet.excluded_calculators = cms.vstring(
	'DileptonCalc',
	'StopCalc',
	'PdfCalc',
	'ChargedHiggsCalc',
	'TprimeCalc',
	'LjetsTopoCalc',
	'LjetsTopoCalcNew',
	'WprimeCalc',
	'JetSubCalc'
	) 

# common calculator options
process.load('LJMet.Com.commonCalc_cfi')

# Stop calculator options
process.load('LJMet.Com.stopCalc_cfi')

# singleLep calculator options
process.load('LJMet.Com.singleLepCalc_cfi')

# LjetsTopoCalc options
#process.load('LJMet.Com.ljetsTopoCalcNew_cfi')
#process.LjetsTopoCalcNew.useBestTop = cms.bool(True)

# Stop calculator options
#process.load('LJMet.Com.JetSubCalc_cfi')


############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('singleLepSelector'),

    # cuts
    debug  = cms.bool(False),

    isMc  = cms.bool(True),
    keepFullMChistory = cms.bool(True),
    
    trigger_cut  = cms.bool(True),
    dump_trigger = cms.bool(False),
    
    mctrigger_path_el = cms.string('HLT_Ele32_eta2p1_WP85_Gsf_v1'),
    mctrigger_path_mu = cms.string('HLT_IsoMu24_eta2p1_IterTrk02_v1'),
    #mctrigger_path_el = cms.string('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1'),
    #mctrigger_path_mu = cms.string('HLT_Mu40_eta2p1_PFJet200_PFJet50_v1'),
    trigger_path_el = cms.vstring('HLT_Ele27_WP80_v8','HLT_Ele27_WP80_v9','HLT_Ele27_WP80_v10','HLT_Ele27_WP80_v11'),
    trigger_path_mu = cms.vstring('HLT_IsoMu24_eta2p1_v11','HLT_IsoMu24_eta2p1_v12','HLT_IsoMu24_eta2p1_v13','HLT_IsoMu24_eta2p1_v14','HLT_IsoMu24_eta2p1_v15'),

    pv_cut         = cms.bool(True),
    hbhe_cut       = cms.bool(True),

    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(2.4),
    min_jet                  = cms.int32(2),
    max_jet                  = cms.int32(4000),
    leading_jet_pt	     = cms.double(100.),

    muon_cuts                = cms.bool(True),
    muon_selector            = cms.bool(False),
    muon_minpt               = cms.double(35.0),
    muon_reliso              = cms.double(0.12),
    muon_maxeta              = cms.double(2.1),
    loose_muon_selector      = cms.bool(False),
    loose_muon_selector_tight= cms.bool(False),
    loose_muon_minpt         = cms.double(25.0),
    loose_muon_reliso        = cms.double(0.2),
    loose_muon_maxeta        = cms.double(2.4),
    min_muon                 = cms.int32(0),

    electron_cuts            = cms.bool(True),
    electron_minpt           = cms.double(35.0),
    electron_maxeta          = cms.double(2.5),
    loose_electron_minpt     = cms.double(25.0),
    loose_electron_maxeta    = cms.double(2.5),
    min_electron             = cms.int32(0),
    
    min_lepton               = cms.int32(1),
    max_lepton               = cms.int32(1),
    second_lepton_veto       = cms.bool(True),
    #trigger_consistent       = cms.bool(False),
    tau_veto	             = cms.bool(False),
    
    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(20.0),
    
    btag_cuts                = cms.bool(False),
    btagger                  = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    btag_min_discr           = cms.double(0.814),
    btag_1                   = cms.bool(False),
    btag_2                   = cms.bool(False),
    btag_3                   = cms.bool(False),

    BTagUncertUp             = cms.bool(False),
    BTagUncertDown           = cms.bool(False),
    JECup                    = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    JEC_txtfile              = cms.string('$CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V5_DATA_UncertaintySources_AK5PF.txt'),
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    jet_collection           = cms.InputTag('slimmedJets'),
    muon_collection          = cms.InputTag('slimmedMuons'),
    electron_collection      = cms.InputTag('slimmedElectrons'),
    tau_collection	     = cms.InputTag('slimmedTaus'),
    met_collection           = cms.InputTag('slimmedMETs'),
#    type1corrmet_collection  = cms.InputTag(''),
    doNewJEC		     = cms.bool(True),
    doLepJetCleaning	     = cms.bool(True),

    MCL1JetPar               = cms.string("../data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt"),
    MCL2JetPar               = cms.string("../data/PHYS14_25_V2_L2Relative_AK4PFchs.txt"),
    MCL3JetPar               = cms.string("../data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt"),

    MCL1JetParAK8            = cms.string("../data/PHYS14_25_V2_L1FastJet_AK8PFchs.txt"),
    MCL2JetParAK8            = cms.string("../data/PHYS14_25_V2_L2Relative_AK8PFchs.txt"),
    MCL3JetParAK8            = cms.string("../data/PHYS14_25_V2_L3Absolute_AK8PFchs.txt"),

    DataL1JetPar             = cms.string('$CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L1FastJet_AK5PFchs.txt'),
    DataL2JetPar             = cms.string('$CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L2Relative_AK5PFchs.txt'),
    DataL3JetPar             = cms.string('$CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L3Absolute_AK5PFchs.txt'),
    DataResJetPar            = cms.string('$CMSSW_BASE/src/LJMet/singletPrime/JEC/Summer13_V4_DATA_L2L3Residual_AK5PFchs.txt')
    )


#######################################################
#
# Input files
#

process.inputs = cms.PSet (
       nEvents    = cms.int32(100000),
           skipEvents = cms.int32(0),
           lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
           fileNames  = cms.vstring(
			'file:///uscms_data/d3/drankin/Wprimetotb/gen/CMSSW_7_2_0/src/PHYS14_Step3/SingletopWprime_M2700.root'
                   )
       )

#process.load(input_module)


# JSON
if (not process.ljmet.isMc==cms.bool(True)):
    JsonFile = 'CMSSW_BASE/src/LJMet/singletPrime/json/Jan222013ReReco_json.txt'
    myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
    process.inputs.lumisToProcess.extend(myList)
        
                
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('Wprime2700Right_singleLep'),
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
process.load('LJMet.Com.pfMuonSelector_cfi') 

# Loose muon
process.LoosepfMuonSelector = process.pfMuonSelector.clone()

# Tight electron
process.load('LJMet.Com.TopElectronSelector_cfi')
#process.TopElectronSelector.cutsToIgnore = cms.vstring("reliso_EB","reliso_EE")

# Loose electron
process.LooseTopElectronSelector = process.TopElectronSelector.clone()
process.LooseTopElectronSelector.version = cms.string('VETO')
