import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

# Define base process
process = cms.Process("LJMetCom")

# Arguments from condor submit script which are used more than once -- EDIT USER DIRECTORY
condorIsMC = bool(True)
relBase    = str('/uscms/home/jmanagan/work/CMSSW_7_3_0/')
condorJSON = str('Jan222013ReReco_json.txt')

# Dilepton calculator options
process.load('LJMet.Com.DileptonCalc_cfi')
process.DileptonCalc.isMc     = condorIsMC
process.DileptonCalc.dataType = cms.string('All')

############################################################
#
# FWLite application options
#
process.ljmet = cms.PSet(
    isMc = cms.bool(condorIsMC),
    runs = cms.vint32([])
)

#Exclude unnecessary calculators
process.ljmet.excluded_calculators = cms.vstring(
    'WprimeCalc',
    'singleLepCalc',
    'LjetsTopoCalc',
    'LjetsTopoCalcNew',
    'StopCalc',
    'PdfCalc',
    'TprimeCalc',
    'ChargedHiggsCalc',
    'TopEventReweightCalc'
    )
############################################################
#
# Dilepton selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('DileptonSelector'),
    isMc      = cms.bool(condorIsMC),

    # Define the cuts -- variable names are strings searched by src/DileptonEventSelector.cc

    # HLT
    trigger_cut              = cms.bool(True),
    dump_trigger             = cms.bool(False),

    # Can use same trigger paths for data and MC since MC is always one of the data versions
    # SOME TRIGGERS ARE LIKELY OLD
    trigger_path_ee          = cms.vstring('HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1', 'HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v1',
                                           'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15',
                                           'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16',
                                           'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17',
                                           'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18',
                                           'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19'),
    
    trigger_path_em          = cms.vstring('HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1', 'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1',
                                           'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
                                           'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
                                           'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8', 'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9',
                                           'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
                                           'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
                                           'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8', 'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9'),
    
    trigger_path_mm          = cms.vstring('HLT_Mu17_Mu8_v1', 'HLT_Mu17_TkMu8_v1', 'HLT_Mu30_TkMu11_v1', 
                                           'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1',
                                           'HLT_Mu17_Mu8_v16', 'HLT_Mu17_Mu8_v17', 'HLT_Mu17_Mu8_v18',
                                           'HLT_Mu17_Mu8_v19', 'HLT_Mu17_Mu8_v21', 'HLT_Mu17_Mu8_v22',
                                           'HLT_Mu17_TkMu8_v9',  'HLT_Mu17_TkMu8_v10', 'HLT_Mu17_TkMu8_v11',
                                           'HLT_Mu17_TkMu8_v12', 'HLT_Mu17_TkMu8_v13', 'HLT_Mu17_TkMu8_v14'),    

    # PV related cuts
    pv_cut                   = cms.bool(True),
    hbhe_cut                 = cms.bool(True),

    # Jet related cuts
    jet_cuts                 = cms.bool(False),
    jet_minpt                = cms.double(20.0),
    jet_maxeta               = cms.double(5),
    min_jet                  = cms.int32(0),
    max_jet                  = cms.int32(4000),

    # Muon related cuts
    muon_cuts                = cms.bool(True),
    min_muon                 = cms.int32(0),
    muon_minpt               = cms.double(10.0),
    muon_maxeta              = cms.double(2.5),
    max_muon                 = cms.int32(20),

    # Electron related cuts
    electron_cuts            = cms.bool(True),
    min_electron             = cms.int32(0),
    electron_minpt           = cms.double(10.0),
    electron_maxeta          = cms.double(2.5),
    max_electron             = cms.int32(20),

    min_lepton               = cms.int32(2),

    # MET related cuts
    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(0.0),

    # Btagging cuts
    btag_cuts                = cms.bool(False),
    btagOP                   = cms.string("CSVL"),
    btag_1                   = cms.bool(True),
    btag_2                   = cms.bool(True),
    btag_3                   = cms.bool(False),

    # Define the branch names of object collections in the input miniAOD file
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    jet_collection           = cms.InputTag('slimmedJets'),
    muon_collection          = cms.InputTag('slimmedMuons'),
    electron_collection      = cms.InputTag('slimmedElectrons'),
    met_collection           = cms.InputTag('slimmedMETs'),

    # Jet corrections are possible, read from txt files. These will certainly need updating!
    JEC_txtfile = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V5_DATA_UncertaintySources_AK5PF.txt'),
    JECup		     = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),

    MCL1JetPar = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V4_MC_L1FastJet_AK5PF.txt'),
    MCL2JetPar = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V4_MC_L2Relative_AK5PF.txt'),
    MCL3JetPar = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V4_MC_L3Absolute_AK5PF.txt'),

    DataL1JetPar = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V4_DATA_L1FastJet_AK5PF.txt'),
    DataL2JetPar = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V4_DATA_L2Relative_AK5PF.txt'),
    DataL3JetPar = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V4_DATA_L3Absolute_AK5PF.txt'),
    DataResJetPar = cms.string(relBase+'/src/LJMet/singleLepton/JEC/Summer13_V4_DATA_L2L3Residual_AK5PF.txt')
)

#######################################################
#
# Input files
#

process.inputs = cms.PSet (
    nEvents    = cms.int32(1000),
    skipEvents = cms.int32(0),
    useHcalLaserEventFilter = cms.bool(True),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),

    # files can be local or accessed from the grid, in which case you need a proxy
    fileNames = cms.vstring('root://cmsxrootd-site.fnal.gov//store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/52DA156A-8470-E411-A8BF-0025905A60AA.root',)
)

# JSON
if not condorIsMC:
    JsonFile = relBase+'/src/LJMet/singleLepton/json/'+condorJSON
    myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
    process.inputs.lumisToProcess.extend(myList)
        
        
        
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('ljmet_dilep'),
    treeName   = cms.string('ljmet'),
)

#######################################################
#
# Object selector options (also see DileptonCalc_cfi.py for the lepton selectors)
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
process.pfJetIDSelector.quality = cms.string('TIGHT')
