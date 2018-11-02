/*
 Contact:        Sadia Khalil (skhalil@fnal.gov)
 */
#ifndef PhysicsTools_PatUtils_interface_TopElectronSelector_h
#define PhysicsTools_PatUtils_interface_TopElectronSelector_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//Math
#include "CLHEP/Units/GlobalPhysicalConstants.h"
typedef math::XYZPoint Point;

class TopElectronSelector : public Selector<pat::Electron>  {
    
public: // interface
    
    bool verbose_;
    bool runData_;
    
    void setUseData(const bool &flag) { runData_ = flag; }
    enum Version_t { VETO, LOOSE, MEDIUM, TIGHT, NONE, HEEP, N_VERSIONS};
    TopElectronSelector() {}
    
    
    TopElectronSelector( edm::ParameterSet const & parameters ){
        
        verbose_ = true;
        
        //runData_= parameters.getParameter<bool>("runData");
        std::string versionStr = parameters.getParameter<std::string>("version");
        Version_t version = N_VERSIONS;
        
        if (versionStr == "VETO"){
            version = VETO;
            if(verbose_) std::cout << "TopElectronSelector: You have choosen version = VETO " <<  std::endl;
        }
        
        if (versionStr == "LOOSE"){
            version = LOOSE;
            if(verbose_) std::cout << "TopElectronSelector: You have choosen version = LOOSE " <<  std::endl;
        }
        
        if (versionStr == "MEDIUM"){
            version = MEDIUM;
            if(verbose_) std::cout << "TopElectronSelector: You have choosen version = MEDIUM " <<  std::endl;
        }
        
        if (versionStr == "TIGHT"){
            version = TIGHT;
            if(verbose_) std::cout << "TopElectronSelector: You have choosen version = TIGHT " <<  std::endl;
        }
        
        if ( versionStr == "NONE" ){
            version = NONE;
            if(verbose_){ std::cout << "TopElectronSelector: If you want to use version NONE "
                << "then make sure to provide the selection cuts by yourself " << std::endl;}
        }

        if (versionStr == "HEEP"){
            version = HEEP;
            if(verbose_) {std::cout << "TopElectronSelector: You have choosen version = HEEP which is non-configurable from ljmet, the cuts are hardcoded in TopElectronSelector. Be Careful!" << std::endl;}
        }
        
        initialize( version,
                   parameters.getParameter<Double_t>("deta_EB"),
                   parameters.getParameter<Double_t>("dphi_EB"),
                   parameters.getParameter<Double_t>("sihih_EB"),
                   parameters.getParameter<Double_t>("hoe_EB"),
                   parameters.getParameter<Double_t>("d0_EB"),
                   parameters.getParameter<Double_t>("dZ_EB"),
                   parameters.getParameter<Double_t>("ooemoop_EB"),
                   parameters.getParameter<Double_t>("reliso_EB"),
                   parameters.getParameter<Double_t>("deta_EE"),
                   parameters.getParameter<Double_t>("dphi_EE"),
                   parameters.getParameter<Double_t>("sihih_EE"),
                   parameters.getParameter<Double_t>("hoe_EE"),
                   parameters.getParameter<Double_t>("d0_EE"),
                   parameters.getParameter<Double_t>("dZ_EE"),
                   parameters.getParameter<Double_t>("ooemoop_EE"),
                   parameters.getParameter<Double_t>("reliso_EE"),
                   parameters.getParameter<Int_t>("mHits_EB"),
                   parameters.getParameter<Int_t>("mHits_EE"),
                   parameters.getParameter<Bool_t>("vtxFitConv")


                   );
        
        
        
        if ( parameters.exists("cutsToIgnore") )
            setIgnoredCuts( parameters.getParameter<std::vector<std::string> >("cutsToIgnore") );
        
        retInternal_ = getBitTemplate();
        pvSrc_  = parameters.getParameter<edm::InputTag>("pvSrc");
        rhoSrc_ = parameters.getParameter<edm::InputTag>("rhoSrc");
        
    }
    
    void initialize(Version_t version,
                    Double_t deta_EB, Double_t  dphi_EB, Double_t sihih_EB, Double_t hoe_EB, Double_t d0_EB, Double_t dZ_EB, Double_t ooemoop_EB, Double_t reliso_EB,
                    Double_t deta_EE, Double_t  dphi_EE, Double_t sihih_EE, Double_t hoe_EE, Double_t d0_EE, Double_t dZ_EE, Double_t ooemoop_EE, Double_t reliso_EE,
                    Int_t mHits_EB, Int_t mHits_EE, Bool_t vtxFitConv)
    {
        version_ = version;
        
        push_back("deta_EB"    );
        push_back("dphi_EB"    );
        push_back("sihih_EB"   );
        push_back("hoe_EB"     );
        push_back("d0_EB"      );
        push_back("dZ_EB"      );
        push_back("ooemoop_EB" );
        push_back("reliso_EB"  );
        push_back("deta_EE"    );
        push_back("dphi_EE"    );
        push_back("sihih_EE"   );
        push_back("hoe_EE"     );
        push_back("d0_EE"      );
        push_back("dZ_EE"      );
        push_back("ooemoop_EE" );
        push_back("reliso_EE"  );
        push_back("mHits_EB"      );
        push_back("mHits_EE"      );
        push_back("vtxFitConv" );
        

        if (version_ == NONE){
            set("deta_EB",     deta_EB);
            set("dphi_EB",     dphi_EB);
            set("sihih_EB",    sihih_EB);
            set("hoe_EB",      hoe_EB);
            set("d0_EB",       d0_EB);
            set("dZ_EB",       dZ_EB);
            set("ooemoop_EB",  ooemoop_EB);
            set("reliso_EB",   reliso_EB);
            set("deta_EE",     deta_EE);
            set("dphi_EE",     dphi_EE);
            set("sihih_EE",    sihih_EE);
            set("hoe_EE",      hoe_EE);
            set("d0_EE",       d0_EE);
            set("dZ_EE",       dZ_EE);
            set("ooemoop_EE",  ooemoop_EE);
            set("reliso_EE",   reliso_EE);
            set("mHits_EB",    mHits_EB);
            set("mHits_EE",    mHits_EE);
            set("vtxFitConv",  vtxFitConv);
        }
        
        if (version_ == VETO) {
            set("deta_EB",     0.013625);
            set("dphi_EB",     0.230374);
            set("sihih_EB",    0.011586);
            set("hoe_EB",      0.181130);
            set("d0_EB",       0.094095);
            set("dZ_EB",       0.713070);
            set("ooemoop_EB",  0.295751);
            set("reliso_EB",   0.158721);
            set("deta_EE",     0.011932);
            set("dphi_EE",     0.255450);
            set("sihih_EE",    0.031849);
            set("hoe_EE",      0.223870);
            set("d0_EE",       0.342293);
            set("dZ_EE",       0.953461);
            set("ooemoop_EE",  0.155501);
            set("reliso_EE",   0.177032);
            set("mHits_EB",    2);
            set("mHits_EE",    3);
            set("vtxFitConv",  0);
        }
        
        if (version_ == LOOSE) {
            set("deta_EB",     0.009277);
            set("dphi_EB",     0.094739);
            set("sihih_EB",    0.010331);
            set("hoe_EB",      0.093068);
            set("d0_EB",       0.035904);
            set("dZ_EB",       0.075496);
            set("ooemoop_EB",  0.189968);
            set("reliso_EB",   0.130136);
            set("deta_EE",     0.009833);
            set("dphi_EE",     0.149934);
            set("sihih_EE",    0.031838);
            set("hoe_EE",      0.115754);
            set("d0_EE",       0.099266);
            set("dZ_EE",       0.197897);
            set("ooemoop_EE",  0.140662);
            set("reliso_EE",   0.163368);
            set("mHits_EB",    1);
            set("mHits_EE",    1);
            set("vtxFitConv",  0);
        }
        
        if (version_ == MEDIUM) {
            set("deta_EB",     0.008925);
            set("dphi_EB",     0.035973);
            set("sihih_EB",    0.009996);
            set("hoe_EB",      0.050537);
            set("d0_EB",       0.012235);
            set("dZ_EB",       0.042020);
            set("ooemoop_EB",  0.091942);
            set("reliso_EB",   0.097213);
            set("deta_EE",     0.007429);
            set("dphi_EE",     0.067879);
            set("sihih_EE",    0.030135);
            set("hoe_EE",      0.086782);
            set("d0_EE",       0.036719);
            set("dZ_EE",       0.138142);
            set("ooemoop_EE",  0.100683);
            set("reliso_EE",   0.113254);
            set("mHits_EB",    1);
            set("mHits_EE",    1);
            set("vtxFitConv",  0);
        }
        

        if (version_ == TIGHT) {
            set("deta_EB",     0.006046);
            set("dphi_EB",     0.028092);
            set("sihih_EB",    0.009947);
            set("hoe_EB",      0.045772);
            set("d0_EB",       0.008790);
            set("dZ_EB",       0.021226);
            set("ooemoop_EB",  0.020118);
            set("reliso_EB",   0.069537);
            set("deta_EE",     0.007057);
            set("dphi_EE",     0.030159);
            set("sihih_EE",    0.028237);
            set("hoe_EE",      0.067778);
            set("d0_EE",       0.027984);
            set("dZ_EE",       0.133431);
            set("ooemoop_EE",  0.098919);
            set("reliso_EE",   0.078265);
            set("mHits_EB",    1);
            set("mHits_EE",    1);
            set("vtxFitConv",  0);
        }
        
        indexDphi_EB_       = index_type(&bits_, "dphi_EB"      );
        indexDeta_EB_       = index_type(&bits_, "deta_EB"      );
        indexSinhih_EB_     = index_type(&bits_, "sihih_EB"     );
        indexHoE_EB_        = index_type(&bits_, "hoe_EB"       );
        indexD0_EB_         = index_type(&bits_, "d0_EB"        );
        indexDZ_EB_         = index_type(&bits_, "dZ_EB"        );
        indexOoemoop_EB_    = index_type(&bits_, "ooemoop_EB"   );
        indexRelIso_EB_     = index_type(&bits_, "reliso_EB"    );
        indexDphi_EE_       = index_type(&bits_, "dphi_EE"      );
        indexDeta_EE_       = index_type(&bits_, "deta_EE"      );
        indexSinhih_EE_     = index_type(&bits_, "sihih_EE"     );
        indexHoE_EE_        = index_type(&bits_, "hoe_EE"       );
        indexD0_EE_         = index_type(&bits_, "d0_EE"        );
        indexDZ_EE_         = index_type(&bits_, "dZ_EE"        );
        indexOoemoop_EE_    = index_type(&bits_, "ooemoop_EE"   );
        indexRelIso_EE_     = index_type(&bits_, "reliso_EE"    );
        indexMHits_EB_      = index_type(&bits_, "mHits_EB"     );
        indexMHits_EE_      = index_type(&bits_, "mHits_EE"     );
        indexVtxFitConv_    = index_type(&bits_, "vtxFitConv"   );
    }
    
    using Selector<pat::Electron>::operator();
    
    
    // Allow for multiple definitions of the cuts.
    
    bool operator()( const pat::Electron & electron, edm::EventBase const & event, pat::strbitset & ret)
    {
        edm::Handle<std::vector<reco::Vertex> > pvtxHandle;
        event.getByLabel( pvSrc_, pvtxHandle);
        Int_t PVsize = pvtxHandle->size();
        if ( PVsize > 0 ) {
            PVtx = pvtxHandle->at(0).position();
        } else {
            throw cms::Exception("InvalidInput") << " There needs to be at least one primary vertex in the event." << std::endl;
        }
        
        edm::Handle<double> rhoHandle;
        event.getByLabel(rhoSrc_, rhoHandle);
        rhoIso = std::max(*(rhoHandle.product()), 0.0);
        
	edm::InputTag convLabel_ ("reducedEgamma:reducedConversions");
	event.getByLabel(convLabel_, conversions);
        edm::InputTag bsLabel_ ("offlineBeamSpot");
	event.getByLabel(bsLabel_, bsHandle);

        return operator()(electron, ret);
    }
    
    //Spring 12 cuts
    bool operator()( const pat::Electron & electron, pat::strbitset & ret)
    {
        ret.set(false);
	if (version_ == HEEP) {
	    Double_t heepPt = electron.ecalDrivenMomentum().pt();
	    Double_t heepEta = electron.superCluster()->eta();
            Int_t heepEcalDriven = electron.ecalDriven();
	    Double_t heepDetainSeed = electron.deltaEtaSeedClusterTrackAtVtx();
	    Double_t heepDphiin = electron.deltaPhiSuperClusterTrackAtVtx();
	    Double_t heepSihih = electron.full5x5_sigmaIetaIeta();
	    Double_t heepE2x5oE5x5 = electron.full5x5_e5x5()!=0 ? electron.full5x5_e2x5Max()/electron.full5x5_e5x5() : 0;
	    Double_t heepE1x5oE5x5 = electron.full5x5_e5x5()!=0 ? electron.full5x5_e1x5()/electron.full5x5_e5x5() : 0;
	    Double_t heepEnergy = electron.superCluster()->energy();
	    Double_t heepHoE = electron.hadronicOverEm();
	    Double_t heepIsoTrkPt = electron.dr03TkSumPt();
	    Double_t heepIsoEmHadDepth1 = electron.dr03EcalRecHitSumEt() + electron.dr03HcalDepth1TowerSumEt();
	    Double_t heepEt = electron.et();
	    Double_t heepDxy = ( PVsize ? electron.gsfTrack()->dxy(PVtx) : electron.gsfTrack()->dxy() );
            Int_t heepMHits   =  electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
	    bool heepPass = false;

            if ( heepPt > 35. && fabs(heepEta) < 1.4442 ) {
		while (1) {
                    if (heepEcalDriven == 1) {}
                    else break;
		    if (fabs(heepDetainSeed) < 0.004) {}
		    else break;
		    if (fabs(heepDphiin) < 0.06) {}
		    else break;
		    if (heepHoE*heepEnergy < 0.05*heepEnergy + 1.) {}
		    else break;
		    if (heepE2x5oE5x5 > 0.94 || heepE1x5oE5x5 > 0.83) {}
		    else break;
		    if (heepIsoEmHadDepth1 < 2. + 0.03*heepEt + 0.28*rhoIso) {}
		    else break;
		    if (heepIsoTrkPt < 5.) {}
		    else break;
		    if (heepMHits <= 1) {}
		    else break;
		    if (fabs(heepDxy) < 0.02) {}
		    else break;
		    heepPass = true;
		    break;
		}
	    }
	    else if ( heepPt > 35. && fabs(heepEta) > 1.566 && fabs(heepEta) < 2.5 ) {
		while (1) {
                    if (heepEcalDriven == 1) {}
                    else break;
		    if (fabs(heepDetainSeed) < 0.006) {}
		    else break;
		    if (fabs(heepDphiin) < 0.06) {}
		    else break;
		    if (heepHoE*heepEnergy < 0.05*heepEnergy + 5.) {}
		    else break;
		    if (heepSihih < 0.03) {}
		    else break;
		    if (heepEt < 50.) {
		        if (heepIsoEmHadDepth1 < 2.5 + 0.28*rhoIso) {}
		        else break;
		    }
		    else {
		        if (heepIsoEmHadDepth1 < 2.5 + 0.03*(heepEt-50) + 0.28*rhoIso) {}
		        else break;
		    }
		    if (heepIsoTrkPt < 5.) {}
		    else break;
		    if (heepMHits <= 1) {}
		    else break;
		    if (fabs(heepDxy) < 0.05) {}
		    else break;
		    heepPass = true;
		    break;
		}
            }

	    ret.set(heepPass);
	}
	else {       
            Double_t scEta = electron.superCluster()->eta();
            Double_t AEff;
            if(fabs(scEta) >2.2) AEff = 0.1337;
            else if(fabs(scEta) >2.0) AEff = 0.0727;
            else if(fabs(scEta) >1.3) AEff = 0.0632;
            else if(fabs(scEta) >0.8) AEff = 0.0954;
            else if(fabs(scEta) >=0.0) AEff = 0.0973;

            
            Double_t Deta  = electron.deltaEtaSuperClusterTrackAtVtx();
            Double_t Dphi  = electron.deltaPhiSuperClusterTrackAtVtx();
            Double_t sihih = electron.full5x5_sigmaIetaIeta();
            Double_t HoE   = electron.hcalOverEcal();
            Double_t D0    = (-1.0)*electron.gsfTrack()->dxy(PVtx);
            Double_t DZ    = electron.gsfTrack()->dz(PVtx);//
            
            Double_t Ooemoop;
            if (electron.ecalEnergy()==0) Ooemoop = 999.;
            else if (!std::isfinite(electron.ecalEnergy())) Ooemoop = 998.;
            else Ooemoop = (1.0/electron.ecalEnergy() - electron.eSuperClusterOverP()/electron.ecalEnergy());
            Int_t mHits   =  electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
            //Bool_t vtxFitConv = electron.passConversionVeto();
            const reco::BeamSpot &beamspot = *bsHandle.product();
            Bool_t vtxFitConv = ConversionTools::hasMatchedConversion(electron, conversions, beamspot.position());
            reco::GsfElectron::PflowIsolationVariables pfIso = electron.pfIsolationVariables();
            Double_t RelIso = ( pfIso.sumChargedHadronPt + std::max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rhoIso*AEff) ) / electron.pt();
    	    bool verbosity =false;
    
    	    if (verbosity) {
    	        std::cout << "\tfabs(Deta) = " << fabs(Deta) << std::endl;
    	        std::cout << "\tfabs(Dphi) = " << fabs(Dphi) << std::endl;
    	        std::cout << "\tsihih = " << sihih << std::endl;
    	        std::cout << "\tHoE = " << HoE << std::endl;
    	        std::cout << "\tfabs(D0) = " << fabs(D0) << std::endl;
    	        std::cout << "\tfabs(DZ) = " << fabs(DZ) << std::endl;
    	        std::cout << "\tfabs(Ooemoop) = " << fabs(Ooemoop) << std::endl;
    	        std::cout << "\tRelIso = " << RelIso << std::endl;
    	        std::cout << "\tmHits = " << mHits << std::endl;
    	        std::cout << "\tvtxFitConv = " << vtxFitConv << std::endl;
                std::cout << "\teta = " << scEta << std::endl;
                std::cout << "\tisEB = " << electron.isEB() << std::endl;
                std::cout << "\tisEE = " << electron.isEE() << std::endl;
    	    }
            
            // now apply the cuts
            if (electron.isEB()) { // BARREL case
                // check the EB cuts
                if ( fabs(Deta)    <  cut(indexDeta_EB_,  double()) || ignoreCut(indexDeta_EB_)  ) passCut(ret, indexDeta_EB_);
                else if (verbosity) std::cout<<"failed Deta"<<std::endl;
                if ( fabs(Dphi)    <  cut(indexDphi_EB_,  double()) || ignoreCut(indexDphi_EB_)  ) passCut(ret, indexDphi_EB_);
                else if (verbosity) std::cout<<"failed Dphi"<<std::endl;
                if ( sihih         <  cut(indexSinhih_EB_,double()) || ignoreCut(indexSinhih_EB_)) passCut(ret, indexSinhih_EB_);
                else if (verbosity) std::cout<<"failed sihih"<<std::endl;
                if ( HoE           <  cut(indexHoE_EB_,   double()) || ignoreCut(indexHoE_EB_)   ) passCut(ret, indexHoE_EB_);
                else if (verbosity) std::cout<<"failed HoE"<<std::endl;
                if ( fabs(D0)      <  cut(indexD0_EB_,    double()) || ignoreCut(indexD0_EB_)    ) passCut(ret, indexD0_EB_);
                else if (verbosity) std::cout<<"failed D0"<<std::endl;
                if ( fabs(DZ)      <  cut(indexDZ_EB_,    double()) || ignoreCut(indexDZ_EB_)    ) passCut(ret, indexDZ_EB_);
                else if (verbosity) std::cout<<"failed DZ"<<std::endl;
                if ( fabs(Ooemoop) <  cut(indexOoemoop_EB_, double()) || ignoreCut(indexOoemoop_EB_) ) passCut(ret, indexOoemoop_EB_);
                else if (verbosity) std::cout<<"failed Ooemoop"<<std::endl;
                if ( RelIso        <  cut(indexRelIso_EB_, double()) || ignoreCut(indexRelIso_EB_) ) passCut(ret, indexRelIso_EB_);
                else if (verbosity) std::cout<<"failed RelIso"<<std::endl;
                if ( mHits         <=  cut(indexMHits_EB_, int()) || ignoreCut(indexMHits_EB_) ) passCut(ret, indexMHits_EB_);
                else if (verbosity) std::cout<<"failed mHits"<<std::endl;
                
                // pass all the EE cuts
                passCut(ret, indexDeta_EE_);
                passCut(ret, indexDphi_EE_);
                passCut(ret, indexSinhih_EE_);
                passCut(ret, indexHoE_EE_);
                passCut(ret, indexD0_EE_);
                passCut(ret, indexDZ_EE_);
                passCut(ret, indexOoemoop_EE_);
                passCut(ret, indexRelIso_EE_);
                passCut(ret, indexMHits_EE_);
            } else {  // ENDCAPS case
                // check the EE cuts
                if ( fabs(Deta)    <  cut(indexDeta_EE_,  double()) || ignoreCut(indexDeta_EE_)  ) passCut(ret, indexDeta_EE_);
                else if (verbosity) std::cout<<"failed Deta"<<std::endl;
                if ( fabs(Dphi)    <  cut(indexDphi_EE_,  double()) || ignoreCut(indexDphi_EE_)  ) passCut(ret, indexDphi_EE_);
                else if (verbosity) std::cout<<"failed Dphi"<<std::endl;
                if ( sihih         <  cut(indexSinhih_EE_,double()) || ignoreCut(indexSinhih_EE_)) passCut(ret, indexSinhih_EE_);
                else if (verbosity) std::cout<<"failed sihih"<<std::endl;
                if ( HoE           <  cut(indexHoE_EE_,   double()) || ignoreCut(indexHoE_EE_)   ) passCut(ret, indexHoE_EE_);
                else if (verbosity) std::cout<<"failed HoE"<<std::endl;
                if ( fabs(D0)            <  cut(indexD0_EE_,    double()) || ignoreCut(indexD0_EE_)    ) passCut(ret, indexD0_EE_);
                else if (verbosity) std::cout<<"failed D0"<<std::endl;
                if ( fabs(DZ)            <  cut(indexDZ_EE_,    double()) || ignoreCut(indexDZ_EE_)    ) passCut(ret, indexDZ_EE_);
                else if (verbosity) std::cout<<"failed DZ"<<std::endl;
                if ( fabs(Ooemoop) <  cut(indexOoemoop_EE_, double()) || ignoreCut(indexOoemoop_EE_) ) passCut(ret, indexOoemoop_EE_);
                else if (verbosity) std::cout<<"failed Ooemoop"<<std::endl;
                if ( RelIso        <  cut(indexRelIso_EE_, double()) || ignoreCut(indexRelIso_EE_) ) passCut(ret, indexRelIso_EE_);
                else if (verbosity) std::cout<<"failed RelIso"<<std::endl;
                if ( mHits         <=  cut(indexMHits_EE_, int()) || ignoreCut(indexMHits_EE_) ) passCut(ret, indexMHits_EE_);
                else if (verbosity) std::cout<<"failed mHits"<<std::endl;
    
                // pass all the EB cuts
                passCut(ret, indexDeta_EB_);
                passCut(ret, indexDphi_EB_);
                passCut(ret, indexSinhih_EB_);
                passCut(ret, indexHoE_EB_);
                passCut(ret, indexD0_EB_);
                passCut(ret, indexDZ_EB_);
                passCut(ret, indexOoemoop_EB_);
                passCut(ret, indexRelIso_EB_);
                passCut(ret, indexMHits_EB_);
                
            }
            if (vtxFitConv     ==  cut(indexVtxFitConv_, bool()) || ignoreCut(indexVtxFitConv_) ) passCut(ret, indexVtxFitConv_);
            
            setIgnored(ret);
	}
	//ret.print(std::cout);
        //std::cout<<std::endl;
        return (bool)ret;
    }
    
private: // member variables
    Version_t version_;
    edm::InputTag pvSrc_;
    Int_t PVsize;
    Point PVtx;
    edm::Handle<reco::ConversionCollection> conversions;
    edm::Handle<reco::BeamSpot> bsHandle;
    edm::InputTag rhoSrc_;
    Double_t rhoIso;
    index_type indexDphi_EB_;
    index_type indexDeta_EB_;
    index_type indexSinhih_EB_;
    index_type indexHoE_EB_;
    index_type indexD0_EB_;
    index_type indexDZ_EB_;
    index_type indexOoemoop_EB_;
    index_type indexRelIso_EB_;
    index_type indexDphi_EE_;
    index_type indexDeta_EE_;
    index_type indexSinhih_EE_;
    index_type indexHoE_EE_;
    index_type indexD0_EE_;
    index_type indexDZ_EE_;
    index_type indexOoemoop_EE_;
    index_type indexRelIso_EE_;
    index_type indexMHits_EB_;
    index_type indexMHits_EE_;
    index_type indexVtxFitConv_;
};

#endif
