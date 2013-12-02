import FWCore.ParameterSet.Config as cms

ljmet = cms.PSet(
    isMc      = cms.bool(True),
    verbosity = cms.int32(0),
    runs                 = cms.vint32([]),
    excluded_calculators = cms.vstring()
    )
