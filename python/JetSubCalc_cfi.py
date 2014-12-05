import FWCore.ParameterSet.Config as cms

# available bDiscrimiant:
# - jetBProbabilityBJetTags
# - jetProbabilityBJetTags
# - trackCountingHighPurBJetTags
# - trackCountingHighEffBJetTags
# - simpleSecondaryVertexHighEffBJetTags
# - simpleSecondaryVertexHighPurBJetTags
# - combinedSecondaryVertexBJetTags (this is the default)
# - combinedInclusiveSecondaryVertexBJetTags

JetSubCalc = cms.PSet(
                      slimmedJetColl     = cms.InputTag("slimmedJets"),
                      slimmedJetsAK8Coll = cms.InputTag("slimmedJetsAK8"),
                      bDiscriminant      = cms.string("combinedSecondaryVertexBJetTags"),
                      tagInfo            = cms.string("CATop")
                      )