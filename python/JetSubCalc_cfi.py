import FWCore.ParameterSet.Config as cms
import os
# available bDiscrimiant:
# - combinedInclusiveSecondaryVertexV2BJetTags
# - combinedMVABJetTags
# - jetBProbabilityBJetTags
# - jetProbabilityBJetTags
# - pfCombinedSecondaryVertexBJetTags
# - simpleSecondaryVertexHighEffBJetTags
# - trackCountingHighPurBJetTags
# - trackCountingHighEffBJetTags
# - simpleSecondaryVertexHighPurBJetTags

relBase = os.environ['CMSSW_BASE']

JetSubCalc = cms.PSet(
                      slimmedJetColl     = cms.InputTag("slimmedJets"),
                      slimmedJetsAK8Coll = cms.InputTag("slimmedJetsAK8"),
                      bDiscriminant      = cms.string("pfDeepCSVJetTags:probb"),
                      genParticles       = cms.InputTag("prunedGenParticles"),
                      tagInfo            = cms.string("caTop"),
                      kappa              = cms.double(0.5),
                      useHTT             = cms.bool(False),
                      killHF             = cms.bool(False),
                      doNewJEC           = cms.bool(False),
                      selectedJetsCA15Coll = cms.InputTag("selectedPatJetsCA15PFCHSNoHF"),
                      useL2L3Mass = cms.bool(False),
                      isMc = cms.bool(True),
                      JECup = cms.bool(False),
                      JECdown = cms.bool(False),
                      JERup = cms.bool(False),
                      JERdown = cms.bool(False),
                      MCL2JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt'),
                      MCL3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt'),
                      MCPTResAK8 = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_MC_PtResolution_AK8PFchs.txt'),
                      MCSF = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_MC_SF_AK4PFchs.txt'),
                      DataL2JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt'),
                      DataL3JetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt'),
                      DataResJetParAK8 = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt'),
                      UncertaintyAK8 = cms.string(relBase+'/src/LJMet/Com/data/Spring16_25nsV6_DATA_Uncertainty_AK8PFchs.txt'),
                      puppiCorrPath = cms.string(relBase+'/src/LJMet/Com/PuppiSoftdropMassCorr/weights/puppiCorr.root'),
                      )

#######################################################################
###  
###   NOTE: You must install HTT to set useHTT = true!!
###   In your CMSSW_X_X_X/src/ directory:
###   git cms-merge-topic gkasieczka:htt-v2-74X
###   scramv1 b -r -j 8
###   
#######################################################################
