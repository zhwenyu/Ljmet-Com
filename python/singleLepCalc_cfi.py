import FWCore.ParameterSet.Config as cms

singleLepCalc = cms.PSet(
                         dataType          = cms.string('All'),
                         elTrigMatchFilters = cms.vstring('hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter'),
                         muTrigMatchFilters = cms.vstring('hltL3fL1sMu16orMu25L1f0L2f16QL3Filtered45e2p1Q'),
                         isMc              = cms.bool(True),
                         pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         genParticles = cms.InputTag("prunedGenParticles"),
			 genJets = cms.InputTag("slimmedGenJets"),
			 packedPFCands = cms.InputTag("packedPFCandidates"),
			 triggerCollection = cms.InputTag("TriggerResults::HLT"),
			 triggerSummary = cms.InputTag("selectedPatTrigger"),
                         keepFullMChistory = cms.bool(True),
                         cleanGenJets = cms.bool(True),
                         keepPDGID    = cms.vuint32(1, 2, 3, 4, 5, 6, 21, 11, 12, 13, 14, 15, 16, 24),
                         keepMomPDGID = cms.vuint32(6, 24),
                         rhoSrc            = cms.InputTag("fixedGridRhoFastjetAll"),
                         UseElMVA = cms.bool(True),
                         )
