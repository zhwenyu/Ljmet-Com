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
#tev left single production inclusive:
#    fileNames = cms.vstring(
#)

#tev right single production inclusive:
#    fileNames = cms.vstring(
#)
#700 left pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-700_left/Private_T53T53ToAll_M-700_left_MINIAOD.root"),
#700 right pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-700_right/Private_T53T53ToAll_M-700_right_MINIAOD.root"),
#900 left pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-900_left/Private_T53T53ToAll_M-900_left_MINIAOD.root"),
#900 right pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-900_right/Private_T53T53ToAll_M-900_right_MINIAOD.root"),
#1000 left pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-1000_left/Private_T53T53ToAll_M-1000_left_MINIAOD.root"),
#1000 right pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-1000_right/Private_T53T53ToAll_M-1000_right_MINIAOD.root"),
#1100 left pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-1100_left/Private_T53T53ToAll_M-1100_left_MINIAOD.root"),
#1100 right pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-1100_right/Private_T53T53ToAll_M-1100_right_MINIAOD.root"),
#1300 left pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-1300_left/Private_T53T53ToAll_M-1300_left_MINIAOD.root"),
#1300 right pair-inclusive:
#        fileNames = cms.vstring("file:/eos/uscms/store/user/lpctlbsm/clint/PHYS14/Inclusive_Decays/PU20/T53T53ToAll_M-1300_right/Private_T53T53ToAll_M-1300_right_MINIAOD.root"),



#WW
#fileNames = cms.vstring(
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/06DCE7A0-F93C-E411-9DEA-002590D0B0CC.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/0E8FABAE-F93C-E411-BE7C-00259073E4CA.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/20F72179-F83C-E411-B119-002590D0B0CA.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/285B065B-F83C-E411-8D53-0025907B4ECA.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/2C3AA376-F83C-E411-A3FB-BCAEC53F6D4E.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/34BC6B92-F83C-E411-B41A-00259074AE8A.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/441F9D43-F83C-E411-A75C-0025907B4E18.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/54EB3A50-F93C-E411-97B0-BCAEC53F6D4E.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/6034FA75-F63C-E411-BF9E-002590D0B0CC.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/6C11DD75-F83C-E411-9A20-485B3989720C.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/7264C812-F93C-E411-AC74-002590D0AFD4.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/7A3CB834-FD3C-E411-BE9F-0025907B4F44.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/9486BC2A-F93C-E411-AA61-20CF305B058C.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/981C84AB-F83C-E411-8CB5-525400CE57D5.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/AA8A2779-F83C-E411-AF24-00259073E4EA.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/CCFB75A5-F93C-E411-BFDC-0025907B4E18.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/D88D8EC2-F83C-E411-848E-E0CB4E19F969.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/DAA10D50-F93C-E411-ADBA-002590D0AFEC.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/E64EAF05-F83C-E411-BE04-002590D0B0CC.root',
        #'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/10000/E833FC89-F93C-E411-9CAB-002590D0B056.root',
#)
#ttbar1
#    fileNames  = cms.vstring(
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/007B37D4-8B70-E411-BC2D-0025905A6066.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/007B37D4-8B70-E411-BC2D-0025905A6066.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3ABFDD92-9470-E411-803B-0025905938AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/6E7B863B-8A70-E411-8641-0025905A48F2.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/CEC84E14-8C70-E411-A16E-0025905A497A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/06843FC5-8370-E411-9B8C-0025905A60AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3C61929C-8770-E411-812F-0025905A612C.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/7824C9F8-9570-E411-AE88-0025905A60CE.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/D47855D6-9970-E411-B31D-0025905A612A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0A867F71-8C70-E411-9CC9-0025905A48D6.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3E7EBC6F-8470-E411-BA33-0025905A60B8.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/808C4B6B-8670-E411-ADF4-0025905A60B8.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/DE2881C7-9A70-E411-9F76-0025905A60CA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0C1D0A70-8870-E411-BAB1-0025905A612C.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/405A56D6-B470-E411-A79D-0026189438D6.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/80DE4DC6-8670-E411-9250-0025905B85B2.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/DEA6A8C7-8670-E411-846C-0025905A607A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0EB35257-8470-E411-A458-0025905B85B2.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/466F87C5-8370-E411-9361-0025905A60AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/8AD759C4-8370-E411-A281-0025905B861C.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/DECE7EF9-8870-E411-8747-0025905A607A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/10D92E71-8470-E411-8DE2-0025905A60B8.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/46E1389E-8770-E411-9B7E-0025905A60AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/8CA5C65E-8470-E411-B07B-0025905964BA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/DEF7F5C4-8670-E411-92B5-0025905A612C.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/123FDAC9-8670-E411-8185-0025905A607A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/52DA156A-8470-E411-A8BF-0025905A60AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/94103378-8670-E411-B716-0025905B85B2.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/E20861C7-8670-E411-AD03-0025905A607A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/18EAB3D4-B470-E411-9F8A-0025905A609E.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/58959E8C-8670-E411-AC45-0025905A607A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/98458B0F-8B70-E411-8627-0025905A48FC.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/E8E53CC5-8370-E411-9A69-0025905A60AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/2E47A1C5-8670-E411-9966-0025905B85D8.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/5CFA4174-8870-E411-8FC2-0025905A48F2.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/AA4EFD56-9870-E411-89C9-0025905A60CE.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/F0F39F5C-8C70-E411-9D24-0025905A48D0.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/2EEAB8C4-8670-E411-B71F-0025905A48F2.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/620CC0F8-9570-E411-85D0-0025905A609E.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/B2833A8F-9570-E411-9F45-0025905A60EE.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/F277FA56-9870-E411-B85E-0025905A612E.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/30CDC258-9570-E411-A0CD-0025905A60EE.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/669576F9-8870-E411-BB67-0025905A607A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/B84549C5-8370-E411-BB66-0025905A60AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/FE2C68C6-8670-E411-A1C9-0025905B85B2.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3204089F-9470-E411-83D3-002590593902.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/6A163270-9570-E411-A627-0025905A6094.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/B8C19E87-8670-E411-BF91-0025905A607A.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/38C85669-8470-E411-987B-0025905A60AA.root',
        #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/6AF33073-8470-E411-93E9-0025905A60AA.root',
       #'file:/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/BE82F88A-8A70-E411-8578-0025905A608A.root',
#)
#ttbar2
#    fileNames  = cms.vstring('/eos/uscms/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/F277FA56-9870-E411-B85E-0025905A612E.root')
#DY
    fileNames  = cms.vstring(
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A733A85-7D6C-E411-8D2B-002481E14D72.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/600D5785-7C6C-E411-B90E-002590DBDFE0.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B6F6C960-7B6C-E411-916C-002481E0DE14.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3008BB28-7D6C-E411-AAC2-002590DB91F0.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/629344EC-7C6C-E411-A19B-0025907DC9B0.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B81F3E5F-796C-E411-9105-002590DB91CE.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0EAD09A8-7C6C-E411-B903-0025901D493E.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/34167B14-7E6C-E411-A113-002590DB92A8.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8618D633-7D6C-E411-AB2C-003048F2FE3E.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C27BA5BA-7D6C-E411-BBA9-002590DB9358.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1E4D0DAE-7C6C-E411-B488-002590DB923C.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3A99E6A9-7B6C-E411-ADB4-00266CFFA6F8.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8E36F058-7C6C-E411-8424-0025901D493E.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C84D5C9B-7C6C-E411-8825-002590DB91CE.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2286DCDB-796C-E411-AAB4-002481E14D72.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5610D8D0-7A6C-E411-B3AA-00237DE0BED6.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/94708D15-7E6C-E411-BA0D-002590DB92A8.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CAD84EE9-7C6C-E411-912C-003048D437A0.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2683B2C5-7C6C-E411-BE0B-002590DB9214.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5EC2A65C-7A6C-E411-94D2-002590DB92A8.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/98175E8A-796C-E411-B612-002590DB923C.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F63F9E51-7D6C-E411-AFD9-002590DB92A8.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/28EF4E6A-7D6C-E411-A54F-0025907DCA38.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5EF8B51F-7C6C-E411-B13F-0025907DC9D6.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A266FB5C-796C-E411-B6EE-0025901D493E.root',
        '/eos/uscms/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FEE3CF68-796C-E411-ABF5-002590DB9214.root',
)
    #WZ
#    fileNames = cms.vstring(
#'file:/eos/uscms/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/484D51C6-2673-E411-8AB0-001E67398412.root',  
#'file:/eos/uscms/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/745E51CD-6473-E411-B67C-002590A80DF0.root',  
#'file:/eos/uscms/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C944BF0-2573-E411-9940-002590200840.root',  
#'file:/eos/uscms/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F2BACD92-2573-E411-9720-002590A83192.root',
#'file:/eos/uscms/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5075EDC3-2573-E411-A168-002590A80E08.root',  
#'file:/eos/uscms/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8033DEC0-2673-E411-AD03-001E67398052.root',  
#'file:/eos/uscms/store/mc/Phys14DR/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/ECED17B7-2573-E411-A237-001E67398011.root',
#)
#WJets
#    fileNames = cms.vstring(
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02215B44-2D70-E411-90A3-0025905A60B8.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3A6426CC-3670-E411-B2FF-00248C65A3EC.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/80314180-3670-E411-A2DF-002618943978.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CCBA9B45-2D70-E411-9A67-003048FFCB84.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0603D444-2D70-E411-AF03-002618943922.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4261BF46-2D70-E411-858A-003048FFD744.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/82CA47D8-2D70-E411-9E97-0025905A6136.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CE88D696-3570-E411-848B-0025905A60E4.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08947C88-3570-E411-974E-002618FDA26D.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/46FD8341-2D70-E411-95B4-0025905A60A6.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8A41A67E-3670-E411-B578-0025905A60A8.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D026D5DC-2D70-E411-9BF8-0025905A60C6.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E8B81D9-2D70-E411-94AB-0025905A4964.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4868A78A-3570-E411-B342-002618943924.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8AE5AC41-2D70-E411-A2B2-0025905A60CE.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D405C046-2D70-E411-B366-003048FFD7C2.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1225D443-2D70-E411-9D85-0025905B85F6.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4A1B9E55-3770-E411-A8DF-00261894393C.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C15AC7F-3670-E411-B918-0025905A60AA.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D4BAE949-2D70-E411-BE41-0025905A60D2.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1631FC7F-3670-E411-A264-002618943811.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4CC5DA9B-3570-E411-B115-00259059642E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C5291E3-2D70-E411-89BC-002618943919.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/D6B6BF43-2D70-E411-92C1-002618943961.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1645BE58-3770-E411-881F-002618943854.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5020EF69-3570-E411-AE7F-0025905A610A.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/90A5DEA6-4870-E411-B256-003048FFD736.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DA4E015B-3370-E411-88A4-0025905B8606.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1CC13D57-3770-E411-836D-0026189438E8.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/506C8643-2D70-E411-9954-0026189438CE.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9232A143-2D70-E411-B94D-002354EF3BDB.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DA4F275D-3770-E411-A519-003048FFCB74.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1EC88B3C-2D70-E411-9DBE-0025905A60CA.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/54F5BF45-2D70-E411-84CF-0025905938B4.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/94AF86D9-2D70-E411-BD48-00259059642E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DC679CD7-2D70-E411-A30A-0025905A6068.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/208E8858-3370-E411-A388-00259059642E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/54F7EE98-3570-E411-A208-0025905A6066.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/964A555A-3770-E411-AFCF-002590593876.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DCE71293-3570-E411-989E-0025905A60D0.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/20EDE349-2D70-E411-BC99-0025905A60D2.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/56B5802E-6077-E411-BEC3-0025905B855E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/96552892-3570-E411-9A5D-0025905A612A.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DE595044-2D70-E411-B665-0025905B85D8.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/224CA68C-3570-E411-9E23-0025905A612A.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/60CF8544-2D70-E411-8AC2-0025905A6088.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/98D04843-2D70-E411-B547-0025905A6122.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DE713DCC-3670-E411-924D-0025905A6090.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24CAAA51-2D70-E411-ABE3-00261894388B.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/644C2B5A-5A70-E411-A9E4-0025905B85D6.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9E178352-2D70-E411-B8AD-0025905A60E0.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E001234C-3A70-E411-B6A2-003048FFCB9E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/267961DA-2D70-E411-8D26-0025905B858E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6CCD7248-2D70-E411-A28C-003048FFCB9E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A2981B44-2D70-E411-B2C8-0025905A6060.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EAC75E43-2D70-E411-AC8E-0025905B8592.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/26BC6B51-2D70-E411-9F32-0025905938A4.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6EDF40D9-2D70-E411-8D2D-0025905A48BA.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/AC4CBED8-2D70-E411-A900-002618FDA237.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/ECAA379B-3570-E411-A87F-0025905A6056.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/28F03240-2D70-E411-81C2-0025905A4964.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7249E0CE-3970-E411-9272-0025905A60B0.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/B883C924-D37B-E411-960D-0026189438EA.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F073FAD6-2D70-E411-802D-0025905A60EE.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2ADBD306-2E70-E411-8DD2-0025905A6070.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7673FF67-3F70-E411-BA27-0025905B8590.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BE2712A6-2E70-E411-B699-0026189438A5.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F0FBFFD6-3670-E411-BFC3-003048FFD752.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3211BB5F-5A70-E411-85EB-002618943832.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/781590D6-2D70-E411-8B87-0025905A6056.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BEE396D4-2D70-E411-A356-00261894395B.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F2F0F2C9-3670-E411-9630-00261894390B.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/340688D7-2D70-E411-904B-002618943807.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/787B8091-3570-E411-9A14-0025905A48EC.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C44B90D6-2D70-E411-AADB-002618943947.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F4304C42-2D70-E411-BE99-0025905B8572.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/362100B1-5070-E411-9A13-0025905A60F8.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/78EBA941-2D70-E411-89B0-0025905A60CE.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4E70B7F-3670-E411-8829-0025905A60A8.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FC5BB390-3570-E411-9D77-003048FFCB74.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3681BFDC-2D70-E411-893B-0025905A60C6.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7A490493-3570-E411-97DB-0025905B8562.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C603C646-2D70-E411-B844-003048FFD744.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FE42B27F-3670-E411-875D-0025905A60AA.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/36B937D7-2D70-E411-81E9-0025905A60CE.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7AE385D9-2D70-E411-B1C4-0025905A4964.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/CCA7E399-3570-E411-9701-0025905A612E.root',
#'file:/eos/uscms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FEF24059-3770-E411-A285-00261894391C.root',
#)
#TTZ
#    fileNames = cms.vstring(
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/188E610D-9871-E411-BABD-002481E15008.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/9C990A71-9771-E411-8F63-002590200A88.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C4B6F654-9B71-E411-8219-001E673972C4.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F060F7F3-B171-E411-AE7A-002590200B38.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/76014B9F-B372-E411-B731-0025B3E0654E.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/BC22410E-9871-E411-A161-001E673971C5.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/EC29E084-5D73-E411-B97B-002590200AE0.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/208787DC-8D71-E411-BD05-001E6739672F.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/5CA9BF41-AA72-E411-A46B-002481E150EE.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8A1F89D8-8D71-E411-BFDD-0025B3E06448.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTZJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/FCF68348-8E71-E411-BC8F-002590A37128.root',
#    )
#TTW
#    fileNames = cms.vstring(
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1C772C5B-DD71-E411-B64E-0025904B1370.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/20A8A674-DD71-E411-A04B-0025901D4B22.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64B14A09-CB71-E411-98A3-0025907FD24A.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F697A0CD-C971-E411-9404-003048D3CB3C.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2010CF0F-CB71-E411-9331-002481E0D678.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/5042A9C5-C971-E411-885B-0025901D42C0.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E2C09A49-CB71-E411-9ACE-002590494C92.root',
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/C2923647-C671-E411-86E0-003048D438FE.root',  
#        'file:/eos/uscms/store/mc/Phys14DR/TTWJets_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/E2C2C43B-C671-E411-A551-002590AC4B76.root',
#        )
#1TeV T53 left
#    fileNames = cms.vstring(
#        'file:/eos/uscms/store/mc/Phys14DR/T53TTo2L2Nu_M-1000_left_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/3CB7457A-0B6A-E411-A274-002590D0B094.root',
#        )
#1TeV T53 right
#    fileNames = cms.vstring(
#        '/eos/uscms/store/mc/Phys14DR/T53TTo2L2Nu_M-1000_right_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64C2805D-CA68-E411-9842-002590DBDFE0.root',
#        '/eos/uscms/store/mc/Phys14DR/T53TTo2L2Nu_M-1000_right_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/7E47F049-CA68-E411-B096-0025907FD430.root',
#        '/eos/uscms/store/mc/Phys14DR/T53TTo2L2Nu_M-1000_right_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/18EE8004-CB68-E411-83A1-003048D43980.root',
#        '/eos/uscms/store/mc/Phys14DR/T53TTo2L2Nu_M-1000_right_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/8CA954F2-CA68-E411-B1DA-003048F0E59E.root',
#        )
)


#Don't need a JSON file for MC
# JSON
#JsonFile = '/data1/speer/tblsm/cmssw/CMSSW_5_3_3/src/LJMet/Com/data/json/Cert_190456-202016_8TeV_PromptReco_Collisions12_JSON_MuonPhys.txt'
#myList   = LumiList.LumiList(filename=JsonFile).getCMSSWString().split(',')
#if not process.ljmet.isMc:
#    process.inputs.lumisToProcess.extend(myList)
        
        
        
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
