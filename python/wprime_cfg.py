import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("LJMetCom")

############################################################
#
# FWLite application options
process.load('LJMet.Com.ljmet_cfi')
process.ljmet.isMc = cms.bool(True)

# common calculator options
process.load('LJMet.Com.commonCalc_cfi')
process.CommonCalc.dummy_parameter = cms.string('Dummy parameter value')

# Stop calculator options
process.load('LJMet.Com.stopCalc_cfi')

# Wprime calculator options
process.load('LJMet.Com.wprimeCalc_cfi')
process.WprimeCalc.isWJets = cms.bool(False)

# LjetsTopoCalc options
process.load('LJMet.Com.ljetsTopoCalcNew_cfi')
process.LjetsTopoCalcNew.useBestTop = cms.bool(True)


############################################################
#
# Event selector options
#
process.event_selector = cms.PSet(

    selection = cms.string('WprimeSelector'),

    # cuts
    debug  = cms.bool(False),

    isMc  = cms.bool(True),
    
    trigger_cut  = cms.bool(True),
    dump_trigger = cms.bool(False),
    
    mctrigger_path_el = cms.string('HLT_Ele27_WP80_v10'), 
    mctrigger_path_mu = cms.string('HLT_IsoMu24_eta2p1_v13'), 
    trigger_path_el = cms.vstring('HLT_Ele27_WP80_v8','HLT_Ele27_WP80_v9','HLT_Ele27_WP80_v10','HLT_Ele27_WP80_v11'), 
    trigger_path_mu = cms.vstring('HLT_IsoMu24_eta2p1_v11','HLT_IsoMu24_eta2p1_v12','HLT_IsoMu24_eta2p1_v13','HLT_IsoMu24_eta2p1_v14','HLT_IsoMu24_eta2p1_v15'),

    #testing muons 
    #mctrigger_path_el = cms.string('HLT_Mu40_v12'), 
    #mctrigger_path_mu = cms.string('HLT_Mu40_v12'), 
    #trigger_path_el = cms.vstring('HLT_Mu40_v12'), 
    #trigger_path_mu = cms.vstring('HLT_Mu40_v12'),
    
    pv_cut         = cms.bool(True),
    hbhe_cut       = cms.bool(True),

    jet_cuts                 = cms.bool(True),
    jet_minpt                = cms.double(30.0),
    jet_maxeta               = cms.double(2.4),
    min_jet                  = cms.int32(2),
    max_jet                  = cms.int32(4000),
    leading_jet_pt           = cms.double(50.0),

    muon_cuts                = cms.bool(True),
    tight_muon_minpt         = cms.double(26.0), # 26 GeV
    tight_muon_maxeta        = cms.double(2.1),
    loose_muon_minpt         = cms.double(10.0),
    loose_muon_maxeta        = cms.double(2.4),
    min_tight_muon           = cms.int32(0),

    electron_cuts                = cms.bool(True),
    tight_electron_minpt         = cms.double(30.0), # 30 GeV
    tight_electron_maxeta        = cms.double(2.5),
    loose_electron_minpt         = cms.double(20.0),
    loose_electron_maxeta        = cms.double(2.5),
    min_tight_electron           = cms.int32(0),
    
    min_tight_lepton         = cms.int32(1),
    max_tight_lepton         = cms.int32(1),
    trigger_consistent       = cms.bool(True),
    second_lepton_veto       = cms.bool(True),
    
    met_cuts                 = cms.bool(True),
    min_met                  = cms.double(20.0),
    
    btag_cuts                = cms.bool(True),
    btagger                  = cms.string('combinedSecondaryVertexBJetTags'),
    btag_min_discr           = cms.double(0.679),
    btag_1                   = cms.bool(True),
    btag_2                   = cms.bool(False),
    btag_3                   = cms.bool(False),

    BTagUncertUp             = cms.bool(False),
    BTagUncertDown           = cms.bool(False),
    JECup                    = cms.bool(False),
    JECdown                  = cms.bool(False),
    JERup                    = cms.bool(False),
    JERdown                  = cms.bool(False),
    JEC_txtfile = cms.string('../cond/Summer12_V2_DATA_AK5PF_UncertaintySources.txt'),
    trigger_collection       = cms.InputTag('TriggerResults::HLT'),
    pv_collection            = cms.InputTag('goodOfflinePrimaryVertices'),
    jet_collection           = cms.InputTag('goodPatJetsPFlow'),
    muon_collection          = cms.InputTag('selectedPatMuonsPFlow'),
    electron_collection      = cms.InputTag('selectedPatElectronsPFlow'),
    met_collection           = cms.InputTag('patMETsPFlow'),
    type1corrmet_collection  = cms.InputTag('pfType1CorrectedMet'),

    )


#######################################################
#
# Input files
#
input_module = 'LJMet.Com.Wprime1900Right_cff'
#input_module = 'LJMet.Com.TT_CT10_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v2_TLBSM_53x_v2_sample_cff'
process.load(input_module)
#process.inputs.nEvents    = cms.int32(1000000)
process.inputs.nEvents    = cms.int32(10000)
process.inputs.skipEvents = cms.int32(0)


# JSON
JsonFile = 'data/json/Cert_190456-202016_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt'
myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
if (not process.ljmet.isMc==cms.bool(True)):
    process.inputs.lumisToProcess.extend(myList)
        
                
#######################################################
#
# Output
#
process.outputs = cms.PSet (
    outputName = cms.string('ljmet_tree'),
    treeName   = cms.string('ljmet'),
)


#######################################################
#
# Object selector options
#

# Primary vertex
process.load('PhysicsTools.SelectorUtils.pvSelector_cfi')
process.pvSelector.pvSrc   = cms.InputTag('goodOfflinePrimaryVertices')
process.pvSelector.minNdof = cms.double(4.0)
process.pvSelector.maxZ    = cms.double(24.0)
process.pvSelector.maxRho  = cms.double(2.0)

# Tight muon
process.load('LJMet.Com.pfMuonSelector_cfi') 
process.pfMuonSelector.version          = cms.string('SPRING11')
process.pfMuonSelector.Chi2             = cms.double(10.0)
process.pfMuonSelector.NHits            = cms.int32(0)
process.pfMuonSelector.NValMuHits       = cms.int32(1)
process.pfMuonSelector.D0               = cms.double(0.2)
process.pfMuonSelector.PFIso            = cms.double(0.12) # 0.12
process.pfMuonSelector.nPixelHits       = cms.int32(1)
process.pfMuonSelector.nMatchedStations = cms.int32(2)
process.pfMuonSelector.nLayersWithMeasurement = cms.int32(6)
process.pfMuonSelector.cutsToIgnore     = cms.vstring('TrackerMuon')

# Loose muon
process.looseMuonSelector = process.pfMuonSelector.clone()
process.looseMuonSelector.PFIso        = cms.double(0.2)
process.looseMuonSelector.nLayersWithMeasurement = cms.int32(0) # not sure why it has to be like this
process.looseMuonSelector.cutsToIgnore = cms.vstring('TrackerMuon',
                                                     'Chi2',
                                                     'NHits',
                                                     'NValMuHits',
                                                     'D0',
                                                     'nPixelHits',
                                                     'nMatchedStations',
                                                     #'nLayersWithMeasurement'
                                                     )

# electron
process.load('LJMet.Com.cutbasedIDSelector_cfi')
process.cutbasedIDSelector.version = cms.string('TIGHT')
process.cutbasedIDSelector.cutsToIgnore     = cms.vstring()

# loose electron
process.looseElectronSelector = process.cutbasedIDSelector.clone()
process.looseElectronSelector.version = cms.string('VETO')
process.looseElectronSelector.cutsToIgnore     = cms.vstring()

# jets
process.load('PhysicsTools.SelectorUtils.pfJetIDSelector_cfi') 
process.pfJetIDSelector.version = cms.string('FIRSTDATA')
process.pfJetIDSelector.quality = cms.string('LOOSE')
