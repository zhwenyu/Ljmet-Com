import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("LJMetCom")

#Common calculator options
process.load('LJMet.Com.commonCalc_cfi')

#Dilepton calculator options
process.load('LJMet.Com.DileptonCalc_cfi')


############################################################
#
# FWLite application options
#
process.ljmet = cms.PSet(
    isMc = cms.bool(True),
    runs = cms.vint32([])
)

process.ljmet.excluded_calculators = cms.vstring(
	'LjetsTopoCalcMinPz',
        'LjetsTopoCalcNew'
	'StopCalc',
	'PdfCalc',
	'ChargedHiggsCalc',
	'TprimeCalc',
	'TpTpCalc',
	'LjetsTopoCalc',
	'WprimeCalc',
        'WprimeBoostedCalc'
	) 


process.DileptonCalc.isMc     = process.ljmet.isMc
process.DileptonCalc.dataType = cms.string('ElMu')

############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('DileptonSelector'),
    isMc              = cms.bool(True),
    keepFullMChistory = cms.bool(True),
    # cuts
    #HLT
    trigger_cut              = cms.bool(False),
    dump_trigger             = cms.bool(True),

    #Can use same trigger paths for data and MC since MC is always one of the data versions
    trigger_path_ee          = cms.vstring('HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1', 'HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v1'),
    
    trigger_path_em          = cms.vstring('HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1', 'HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1'),

    trigger_path_mm          = cms.vstring('HLT_Mu17_Mu8_v1', 'HLT_Mu17_TkMu8_v1', 'HLT_Mu30_TkMu11_v1', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1'),    

    pv_cut                   = cms.bool(True),
    hbhe_cut                 = cms.bool(True),

    jet_cuts                 = cms.bool(False),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(2.4),
    min_jet                  = cms.int32(4),
    max_jet                  = cms.int32(4000),

    muon_cuts                = cms.bool(True),
    min_muon                 = cms.int32(0),
    muon_minpt               = cms.double(10.0),
    muon_maxeta              = cms.double(2.5),
    max_muon                 = cms.int32(10),

    electron_cuts            = cms.bool(True),
    min_electron             = cms.int32(0),
    electron_minpt           = cms.double(10.0),
    electron_maxeta          = cms.double(2.5),
    max_electron             = cms.int32(10),

    min_lepton               = cms.int32(2),

    met_cuts                 = cms.bool(False),
    min_met                  = cms.double(0.0),
    btag_cuts                = cms.bool(False),
     btagOP                  = cms.string("CSVM"),
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
    #We will not use it in this exercise
    JEC_txtfile = cms.string('/uscms_data/d2/avetisya/CMSDAS/2015/CMSSW_7_3_0/src/LJMet/Com/tools/Summer13_V5_DATA_Uncertainty_AK5PF.txt'),
    JECup		     = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),

    #new jet energy corrections
    doNewJEC                 = cms.bool(True),

    #triggerstudy info
    doTriggerStudy           = cms.bool(True),
    TriggerBits              = cms.InputTag("TriggerResults","","HLT"),
    TriggerObjects           = cms.InputTag("selectedPatTrigger"),


)

#######################################################
#
# Input files
#


process.inputs = cms.PSet (
    nEvents    = cms.int32(100),
    skipEvents = cms.int32(0),
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    
)


#Don't need a JSON file for MC
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
#    outputName = cms.string('ljmet_tree_X53X53ToAll_M-700_right'),
    outputName = cms.string('ljmet_tree_DYJets'),
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
