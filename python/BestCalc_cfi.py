import FWCore.ParameterSet.Config as cms

BestCalc = cms.PSet(

    #dnnFile = cms.string(relBase+'/src/LJMet/Com/data/BEST_mlp.json'), #
    #    dnnFile = cms.string(relBase+'/src/LJMet/Com/data/BEST_mlp.json'), #                                                                                                                             
    dnnFile = cms.string('/uscms_data/d3/saj32265/CMSSW_9_4_6_patch1/src/LJMet/Com/data/BEST_mlp.json'), #  

    numSubjetsMin = cms.int32(2),
    numDaughtersMin = cms.int32(2),
    jetSoftDropMassMin = cms.double(1.0),
    jetPtMin = cms.double(200.0),
    radiusSmall = cms.double(0.4),
    radiusLarge = cms.double(0.8),
    reclusterJetPtMin = cms.double(20.0),
    jetChargeKappa = cms.double(0.6),
    maxJetSize = cms.int32(4),
    )

