import FWCore.ParameterSet.Config as cms

XConeCalc = cms.PSet(
                      packedPFCandColl_it	 = cms.InputTag("packedPFCandidates"),
                      XConeR			 = cms.double(0.4),
                      XConeNumJets		 = cms.int32(6),
                      XConeBeta			 = cms.double(2.0),
                      usePFchs			 = cms.bool(True),
                      DEBUG			     = cms.bool(False),
                      )
