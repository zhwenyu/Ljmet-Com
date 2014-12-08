import FWCore.ParameterSet.Config as cms

WprimeCalc = cms.PSet(
                      triggerSummary = cms.InputTag("hltTriggerSummaryAOD"),
                      triggerCollection = cms.InputTag("TriggerResults::HLT"),
                      rhoSrc     = cms.InputTag("fixedGridRhoAll", ''),
                      isWJets     = cms.bool(False),
                      isTB       = cms.bool(False)
                      )
