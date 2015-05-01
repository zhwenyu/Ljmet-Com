import FWCore.ParameterSet.Config as cms

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

JetSubCalc = cms.PSet(
                      slimmedJetColl     = cms.InputTag("slimmedJets"),
                      slimmedJetsAK8Coll = cms.InputTag("slimmedJetsAK8"),
                      bDiscriminant      = cms.string("pfCombinedSecondaryVertexBJetTags"),
                      tagInfo            = cms.string("caTop")
                      )