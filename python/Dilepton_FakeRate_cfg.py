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
    isMc = cms.bool(False),
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
    trigger_cut              = cms.bool(True),
    dump_trigger             = cms.bool(True),

    #Can use same trigger paths for data and MC since MC is always one of the data versions
#    trigger_path_ee          = cms.vstring('HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1', 'HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v1'),
    
    trigger_path_e          = cms.vstring('HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v1'),

    trigger_path_m          = cms.vstring('HLT_Mu17_v1'),

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

    min_lepton               = cms.int32(1),

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
    fileNames = cms.vstring(
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/096/00000/22D22D7F-5626-E511-BDE3-02163E011FAB.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/161/00000/7019DC27-9C26-E511-84FF-02163E011CC2.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/162/00000/9CC606D8-4127-E511-8F35-02163E013830.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/163/00000/9C435096-9F26-E511-A1D7-02163E012AB6.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/164/00000/4633CC68-A326-E511-95D0-02163E0124EA.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/167/00000/E661FBF2-A726-E511-8B8B-02163E01207C.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/168/00000/FCB6CB61-BC26-E511-8858-02163E01375B.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/244/00000/084C9A66-9227-E511-91E0-02163E0133F0.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/244/00000/12EE24E2-8F27-E511-80D1-02163E013793.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/251/00000/3A0E6309-8827-E511-B96D-02163E013432.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/252/00000/36D5A4A8-A227-E511-9AAA-02163E0135F3.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/252/00000/6258DF96-A227-E511-BE8A-02163E01259F.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/491/00000/4A5A5D00-C628-E511-ACC7-02163E01414A.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/493/00000/0C69AF3D-CF28-E511-B56A-02163E01414A.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/496/00000/3AA75CED-932C-E511-8248-02163E0133F2.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/497/00000/607EA0EA-E028-E511-BD54-02163E0133FF.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/498/00000/8064CCF6-ED28-E511-87D2-02163E014437.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/499/00000/70310B47-F728-E511-B2EF-02163E0118A2.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/500/00000/0273C876-0E29-E511-8B38-02163E012712.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/521/00000/D28AB765-6629-E511-8AAD-02163E0133D1.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/522/00000/805EB9CD-6129-E511-BF1C-02163E0129A3.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/548/00000/B6D08898-232A-E511-A833-02163E011DDE.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/559/00000/AA62F6DD-A62C-E511-A8EC-02163E013791.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/560/00000/BA599BB8-E129-E511-B26A-02163E0134CC.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/561/00000/5ACDA1DE-FB29-E511-8D8C-02163E0133B5.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/561/00000/CA80E14E-1E2A-E511-8C7D-02163E0122C2.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/30DDF910-5E2A-E511-9F4D-02163E01206A.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/3CE07240-742A-E511-BA88-02163E01258B.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/5CF006D1-602A-E511-95CE-02163E0126E1.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/8EE9BBAA-7E2A-E511-AEF7-02163E0143C0.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/B41B8802-672A-E511-A9EA-02163E012787.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/DCC900B5-972A-E511-9785-02163E012283.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/F0A7C9F3-6B2A-E511-A73B-02163E0126A0.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/562/00000/FE5AD795-6E2A-E511-9C40-02163E012787.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/604/00000/AE22AF42-902A-E511-8A22-02163E012B30.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/612/00000/50DA7894-932A-E511-801E-02163E0136A2.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/628/00000/40EF63A0-B52A-E511-8B57-02163E0133F0.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/638/00000/0CDDB666-E72A-E511-9BFD-02163E011DE5.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/638/00000/B2FC1038-372B-E511-AA94-02163E013481.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/642/00000/2C622272-D02A-E511-9F20-02163E013645.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/0E4B7E28-8D2C-E511-BFDA-02163E01477B.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/1247CF12-932C-E511-B9ED-02163E01354D.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/3A437BCB-912C-E511-96D0-02163E012934.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/7077210E-8F2C-E511-97D5-02163E0138EC.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/9EFCB7EB-C12C-E511-B8BB-02163E012158.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/A42F5B12-9C2C-E511-AAB3-02163E0134C3.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/C2E62796-942C-E511-8869-02163E0121A1.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/CCA6600A-912C-E511-B1EF-02163E0133D0.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/643/00000/D84DA8FC-8F2C-E511-9B5A-02163E01420D.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/716/00000/02C46BE4-302C-E511-A01C-02163E0128BF.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/721/00000/9663EE89-022C-E511-985A-02163E01359E.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/721/00000/CADC920F-E02B-E511-BB9B-02163E01412F.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/781/00000/662647FF-9B2C-E511-85C5-02163E012965.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/00CD59FD-2B2D-E511-8DB2-02163E01267F.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/500D754A-292D-E511-AEBA-02163E0118F2.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/CA2F0F5F-242D-E511-96AB-02163E011D88.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/CAA3B776-312D-E511-923A-02163E01280D.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/883/00000/D030EA86-9C2D-E511-8FCB-02163E014181.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/102/00000/06E11C1F-DA2F-E511-B04A-02163E0117FF.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/116/00000/2EE8BBD0-7730-E511-95F5-02163E013560.root
/store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/252/126/00000/AA0C0594-CA30-E511-BD15-02163E0123C0.root
)
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
