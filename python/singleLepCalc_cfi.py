import FWCore.ParameterSet.Config as cms

singleLepCalc = cms.PSet(
                         dataType          = cms.string('All'),
                         isMc              = cms.bool(True),
                         pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         genParticles = cms.InputTag("prunedGenParticles"),
			 triggerCollection = cms.InputTag("TriggerResults::HLT"),
			 triggerSummary = cms.InputTag("selectedPatTrigger"),
                         keepFullMChistory = cms.bool(True),
                         keepPDGID    = cms.vuint32(1, 2, 3, 4, 5, 21, 11, 12, 13, 14, 15, 16),
                         keepMomPDGID = cms.vuint32(6, 24),
                         rhoSrc            = cms.InputTag("fixedGridRhoAll"),
                         isWJets           = cms.bool(False),
                         isTB              = cms.bool(False),
                         )
