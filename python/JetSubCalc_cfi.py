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
                      killHF             = cms.bool(False),
                      doNewJEC           = cms.bool(False),
                      selectedJetsCA15Coll = cms.InputTag("selectedPatJetsCA15PFCHSNoHF"),
                      useL2L3Mass = cms.bool(False),
                      isMc = cms.bool(True),
                      MCL2JetParAK8 = cms.string('/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV2_MC_L2Relative_AK8PFchs.txt'),
                      MCL3JetParAK8 = cms.string('/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV2_MC_L3Absolute_AK8PFchs.txt'),
                      DataL2JetParAK8 = cms.string('/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_DATA_L2Relative_AK8PFchs.txt'),
                      DataL3JetParAK8 = cms.string('/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_DATA_L3Absolute_AK8PFchs.txt'),
                      DataResJetParAK8 = cms.string('/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_DATA_L2L3Residual_AK8PFchs.txt')
                      )

#######################################################################
###  
###   NOTE: You must install HTT to set useHTT = true!!
###   In your CMSSW_X_X_X/src/ directory:
###   git cms-merge-topic gkasieczka:htt-v2-74X
###   scramv1 b -r -j 8
###   
#######################################################################
