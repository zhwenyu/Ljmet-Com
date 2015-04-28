import FWCore.ParameterSet.Config as cms

TopElectronSelector = cms.PSet(version            = cms.string("TIGHT"),
                   	       deta_EB		  = cms.double(0.1),
                   	       dphi_EB		  = cms.double(0.1),
                   	       sihih_EB		  = cms.double(0.1),
                   	       hoe_EB		  = cms.double(0.1),
                   	       d0_EB		  = cms.double(0.1),
                   	       dZ_EB		  = cms.double(0.1),
                   	       ooemoop_EB	  = cms.double(0.1),
                   	       reliso_EB	  = cms.double(0.1),
                   	       deta_EE		  = cms.double(0.1),
                   	       dphi_EE		  = cms.double(0.1),
                   	       sihih_EE		  = cms.double(0.1),
                   	       hoe_EE		  = cms.double(0.1),
                   	       d0_EE		  = cms.double(0.1),
                   	       dZ_EE		  = cms.double(0.1),
                   	       ooemoop_EE	  = cms.double(0.1),
                   	       reliso_EE	  = cms.double(0.1),
			       mHits_EB		  = cms.int32(1),
			       mHits_EE		  = cms.int32(1),
			       vtxFitConv	  = cms.bool(True),
        		       pvSrc  		  = cms.InputTag("offlineSlimmedPrimaryVertices"),
        		       rhoSrc 		  = cms.InputTag("fixedGridRhoAll"),
                               #cutsToIgnore       = cms.vstring("reliso_EB","reliso_EE")
                               cutsToIgnore       = cms.vstring()
)
