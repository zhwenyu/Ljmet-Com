import FWCore.ParameterSet.Config as cms

MVAElectronSelector = cms.PSet(version            = cms.string("MEDIUM"),
                               cutValue_Bin	  = cms.double(1.),
                               cutValue_Bout	  = cms.double(1.),
                               cutValue_E	  = cms.double(1.),
        		       pvSrc  		  = cms.InputTag("offlineSlimmedPrimaryVertices"),
        		       rhoSrc 		  = cms.InputTag("fixedGridRhoFastjetAll"),
)
