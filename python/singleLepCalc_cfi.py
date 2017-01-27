import FWCore.ParameterSet.Config as cms

singleLepCalc = cms.PSet(
                         dataType          = cms.string('All'),
                         elTrigMatchFilters = cms.vstring('hltEle15VVVLGsfTrackIsoFilter','hltEle32WPTightGsfTrackIsoFilter','hltEle30WPTightGsfTrackIsoFilter'),
                         muTrigMatchFilters = cms.vstring('hltL3MuVVVLIsoFIlter','hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09','hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09','hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q','hltL3fL1sMu25f0TkFiltered50Q'),
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
                         keepPDGIDForce  = cms.vuint32(6,6),
                         keepStatusForce = cms.vuint32(62,22),
                         rhoSrc            = cms.InputTag("fixedGridRhoFastjetAll"),
                         doElSCMETCorr = cms.bool(False),
                         electronCollection = cms.InputTag("slimmedElectrons"),
                         UseElMVA = cms.bool(True),
                         OverrideLHEWeights = cms.bool(False),
                         basePDFname = cms.string("cteq6"),
                         newPDFname = cms.string("PDF4LHC15_nlo_mc_pdfas"),
                         saveLooseLeps = cms.bool(False)
                         )
