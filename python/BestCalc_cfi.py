import FWCore.ParameterSet.Config as cms

BestCalc = cms.PSet(

    dnnFile = cms.string(''),

    numSubjetsMin = cms.int32(2),
    numDaughtersMin = cms.int32(2),
    jetSoftDropMassMin = cms.double(5.0),
    jetPtMin = cms.double(170.0),
    radiusSmall = cms.double(0.4),
    radiusLarge = cms.double(0.8),
    reclusterJetPtMin = cms.double(20.0),
    jetChargeKappa = cms.double(0.6),
    maxJetSize = cms.int32(4),
    )

