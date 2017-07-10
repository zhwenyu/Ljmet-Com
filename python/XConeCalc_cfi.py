import FWCore.ParameterSet.Config as cms

XConeCalc = cms.PSet(
                      packedPFCandColl_it	 = cms.InputTag("packedPFCandidates"),
                      packedGenParticleColl_it	 = cms.InputTag("packedGenParticles"),
                      XConeR			 = cms.double(0.4),
                      XConeNumJets		 = cms.int32(6),
                      VarNumJets		 = cms.bool(True),
                      XConeBeta			 = cms.double(2.0),
                      usePFchs			 = cms.bool(True),
                      doPUPPI			 = cms.bool(True),
                      doGenXCone		 = cms.bool(True),
                      IsMc				 = cms.bool(False),
                      DEBUG			     = cms.bool(False),
                      saveJetConst		 = cms.bool(False),
                      )
