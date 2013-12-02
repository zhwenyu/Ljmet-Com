import FWCore.ParameterSet.Config as cms

from LJMet.Com.cutbasedIDSelector_cfi import *

DileptonCalc = cms.PSet (
    isMc         = cms.bool(False),
    dataType     = cms.string('ElMu'),
    rhoSrc       = cms.InputTag("kt6PFJetsForIsolation", "rho"),
    pvCollection = cms.InputTag("goodOfflinePrimaryVertices"),
    genParticles = cms.InputTag("prunedGenParticles"),
    keepPDGID    = cms.vuint32(8000001, 11, 12, 13, 14, 15, 16),
    keepMomPDGID = cms.vuint32(6, 24, 8000001),
    cutbasedIDSelectorLoose  = cutbasedIDSelector.clone(),
    cutbasedIDSelectorMedium = cutbasedIDSelector.clone(),
    cutbasedIDSelectorTight  = cutbasedIDSelector.clone()
    
)
DileptonCalc.cutbasedIDSelectorLoose.version = cms.string('LOOSE')
DileptonCalc.cutbasedIDSelectorMedium.version = cms.string('MEDIUM')
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
