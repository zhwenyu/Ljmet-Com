import FWCore.ParameterSet.Config as cms

from LJMet.Com.cutbasedIDSelector_cfi import *

DileptonCalc = cms.PSet (
                         isMc         = cms.bool(False),
                         dataType     = cms.string('ElMu'),
#                         rhoSrc       = cms.InputTag("fixedGridRhoAll", "rho"),
                         pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         genParticles = cms.InputTag("prunedGenParticles"),
                         keepPDGID    = cms.vuint32(8000003, 11, 12, 13, 14, 15, 16, 6, 24),
                         keepMomPDGID = cms.vuint32(6, 24, 8000003),
                         cutbasedIDSelectorLoose  = cutbasedIDSelector.clone(),
                         cutbasedIDSelectorMedium = cutbasedIDSelector.clone(),
                         cutbasedIDSelectorTight  = cutbasedIDSelector.clone(),
                         electronsMiniAOD    = cms.InputTag("slimmedElectrons"),                                               
                         eleMVALooseIDMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                         eleMVATightIDMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                         mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                         
                         )
DileptonCalc.cutbasedIDSelectorLoose.version = cms.string('LOOSE')
DileptonCalc.cutbasedIDSelectorMedium.version = cms.string('MEDIUM')
DileptonCalc.cutbasedIDSelectorLoose.mHits_EE = cms.int32(1)
DileptonCalc.cutbasedIDSelectorLoose.mHits_EB = cms.int32(1)
DileptonCalc.cutbasedIDSelectorLoose.cutsToIgnore = cms.vstring('reliso_EB',
                                                                'reliso_EE',
                                                                'd0_EB',
                                                                'd0_EE')
DileptonCalc.cutbasedIDSelectorMedium.cutsToIgnore = cms.vstring('reliso_EB',
                                                                 'reliso_EE',
                                                                 'd0_EB',
                                                                 'd0_EE')
DileptonCalc.cutbasedIDSelectorTight.cutsToIgnore = cms.vstring('reliso_EB',
                                                                'reliso_EE',
                                                                'd0_EB',
                                                                'd0_EE')
