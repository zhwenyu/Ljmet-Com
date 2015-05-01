import FWCore.ParameterSet.Config as cms

CATopoCalc = cms.PSet(
                      debug             = cms.bool(False),
                      AK8slimmedJetColl = cms.InputTag("slimmedJetsAK8")
                      )
