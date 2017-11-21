import FWCore.ParameterSet.Config as cms
import os

relBase    = os.environ['CMSSW_BASE']

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
                      MCL1JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt'),
                      MCL2JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt'),
                      MCL3JetPar               = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt'),
                      DataL1JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt'),
                      DataL2JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt'),
                      DataL3JetPar             = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt'),
                      DataResJetPar            = cms.string(relBase+'/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt'),

                      )
