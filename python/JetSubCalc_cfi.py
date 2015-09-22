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
                      bDiscriminant      = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                      tagInfo            = cms.string("caTop"),
                      kappa              = cms.double(0.5),
                      useHTT             = cms.bool(False),
                      selectedJetsCA15Coll = cms.InputTag("selectedPatJetsCA15PFCHSNoHF")
                      )
