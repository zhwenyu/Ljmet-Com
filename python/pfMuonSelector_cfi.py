import FWCore.ParameterSet.Config as cms

# List of available parameters and their type:
# version, s
# Chi2, d
# minTrackerLayers, i
# minValidMuHits, i
# maxIp, d
# minPixelHits, i
# minMatchedStations, i
# maxZImpact, d
# maxPfRelIso, d
# cutsToIgnore, s

pfMuonSelector = cms.PSet(version            = cms.string('TOPPAG12_LJETS'),
                          Chi2               = cms.double(10.0),
                          minTrackerLayers   = cms.int32(6),
                          minValidMuHits     = cms.int32(1),
                          maxIp              = cms.double(0.2),
                          minPixelHits       = cms.int32(1),
                          minMatchedStations = cms.int32(2),
                          maxZImpact         = cms.double(999999.),
                          maxPfRelIso        = cms.double(0.12),
                          cutsToIgnore       = cms.vstring()
                          )
