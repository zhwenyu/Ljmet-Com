import FWCore.ParameterSet.Config as cms

singleLepCalc = cms.PSet(
                         triggerSummary    = cms.InputTag("hltTriggerSummaryAOD"),
                         triggerCollection = cms.InputTag("TriggerResults::HLT"),
                         rhoSrc            = cms.InputTag("fixedGridRhoAll"),
#                         muonSrc           = cms.InputTag("slimmedMuons"),
                         isWJets           = cms.bool(False),
                         isTB              = cms.bool(False)
                         )
