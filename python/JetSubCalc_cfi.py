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
                      kappa              = cms.double(0.5),
                      killHF             = cms.bool(False),
                      doNewJEC           = cms.bool(False),
                      useL2L3Mass = cms.bool(False),
                      isMc = cms.bool(True),
                      puppiCorrPath = cms.string(relBase+'/src/LJMet/Com/PuppiSoftdropMassCorr/weights/puppiCorr.root'),
                      )
