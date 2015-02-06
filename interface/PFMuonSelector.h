#ifndef PhysicsTools_PatUtils_interface_PFMuonSelector_h
#define PhysicsTools_PatUtils_interface_PFMuonSelector_h

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>

class PFMuonSelector : public Selector<pat::Muon> {
    
public: // interface
    
    bool verbose_;
    
    enum Version_t { TOPPAG12_LJETS, TIGHT_RUN1, TOPPAG12_LJETS_MOD, TOPPAG12_LJETS_VETO, N_VERSIONS };
    
    PFMuonSelector() {}
    
    PFMuonSelector( edm::ParameterSet const & parameters ) {
        
        verbose_ = false;
        
        std::string versionStr = parameters.getParameter<std::string>("version");
        
        Version_t version = N_VERSIONS;
        
        if ( versionStr == "TOPPAG12_LJETS" ) {
            version = TOPPAG12_LJETS;
        }
        else if ( versionStr == "TIGHT_RUN1" ) {
            version = TIGHT_RUN1;
        }
        else if ( versionStr == "TOPPAG12_LJETS_VETO" ) {
            version = TOPPAG12_LJETS_VETO;
        }
        else if ( versionStr == "TOPPAG12_LJETS_MOD" ) {
            version = TOPPAG12_LJETS_MOD;
        }
        else {
            throw cms::Exception("InvalidInput") << "Expect version to be one of TOPPAG12_LJETS, TOPPAG12_LJETS_VETO, TOPPAG12_LJETS_MOD" << std::endl;
        }
        
        initialize( version,
                   parameters.getParameter<double>("Chi2"),
                   parameters.getParameter<int>   ("minTrackerLayers"),
                   parameters.getParameter<int>   ("minValidMuHits"),
                   parameters.getParameter<double>("maxIp"),
                   parameters.getParameter<int>   ("minPixelHits"),
                   parameters.getParameter<int>   ("minMatchedStations"),
                   parameters.getParameter<double>("maxZImpact"),
                   parameters.getParameter<double>("maxPfRelIso")
                   );
        if ( parameters.exists("cutsToIgnore") )
            setIgnoredCuts( parameters.getParameter<std::vector<std::string> >("cutsToIgnore") );
        
        retInternal_ = getBitTemplate();
        
    }
    
    
    
    void initialize( Version_t version,
                    double    chi2             = 10.0,
                    int       minTrackerLayers = 6,
                    int       minValidMuonHits = 1,
                    double    maxIp            = 0.2,
                    int       minPixelHits     = 1,
                    int       minNMatches      = 2,
                    double    maxZImpact       = 0.5,
                    double    pfiso            = 0.12
                    )
    {
        version_ = version;
        
        push_back("GlobalMuon");
        push_back("TrackerMuon");
        push_back("GlobalOrTrackerMuon");
        push_back("Chi2",                chi2   );
        push_back("minTrackerLayers",    minTrackerLayers);
        push_back("minValidMuHits",      minValidMuonHits  );
        push_back("maxIp",               maxIp );
        push_back("maxZImpact",          maxZImpact );
        push_back("minPixelHits",        minPixelHits);
        push_back("minMatchedStations",  minNMatches);
        push_back("maxPfRelIso",         pfiso );
        
        set("GlobalMuon");
        set("TrackerMuon");
        set("GlobalOrTrackerMuon");
        set("Chi2");
        set("minTrackerLayers");
        set("minValidMuHits");
        set("maxIp");
        set("maxZImpact");
        set("minPixelHits");
        set("minMatchedStations");
        set("maxPfRelIso");
        
        indexChi2_             = index_type(&bits_, "Chi2"            );
        indexMinTrackerLayers_ = index_type(&bits_, "minTrackerLayers" );
        indexminValidMuHits_   = index_type(&bits_, "minValidMuHits"      );
        indexMaxIp_            = index_type(&bits_, "maxIp"      );
        indexMaxZImpact_       = index_type(&bits_, "maxZImpact"      );
        indexPixHits_          = index_type(&bits_, "minPixelHits"      );
        indexStations_         = index_type(&bits_, "minMatchedStations");
        indexmaxPfRelIso_      = index_type(&bits_, "maxPfRelIso"           );
        
        
        if (version_ == TOPPAG12_LJETS ){
            set("GlobalMuon",true);
            set("GlobalOrTrackerMuon",false);
            set("TrackerMuon", false);
        }
        
        if (version_ == TIGHT_RUN1 ){ //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon
            set("GlobalMuon",true);
            set("GlobalOrTrackerMuon",false);
            set("TrackerMuon", false);
            
            set("maxPfRelIso", 0.2);
            set("Chi2", 10);
            set("minTrackerLayers", 5);
            set("minValidMuHits", 0);
            set("maxIp", 0.2);
            set("maxZImpact", 0.5);
            set("minPixelHits", 0);
            set("minMatchedStations", 1);
        }
        
        if (version_ == TOPPAG12_LJETS_MOD ){
            set("GlobalMuon",true);
            set("GlobalOrTrackerMuon",false);
            set("TrackerMuon",false);
            set("cutsToIgnore","maxPfRelIso");
        }
        
        if (version_ == TOPPAG12_LJETS_VETO ){
            set("TrackerMuon", false);
            set("GlobalMuon", false);
            
            set("GlobalOrTrackerMuon", true);
            set("maxPfRelIso", 0.2);
            set("Chi2", false);
            set("minTrackerLayers", false);
            set("minValidMuHits", false);
            set("maxIp", false);
            set("maxZImpact", false);
            set("minPixelHits", false);
            set("minMatchedStations", false);
        }
        
    }
    
    // Allow for multiple definitions of the cuts.
    bool operator()( const pat::Muon & muon, pat::strbitset & ret )
    {
        if (version_ == TOPPAG12_LJETS || version_ == TIGHT_RUN1 || version_ == TOPPAG12_LJETS_VETO || version_ == TOPPAG12_LJETS_MOD) return TopPag12LjetsCuts(muon, ret);
        else {
            return false;
        }
    }
    
    
    
    using Selector<pat::Muon>::operator();
    
    
    
    bool TopPag12LjetsCuts( const pat::Muon & muon, pat::strbitset & ret){
        
        ret.set(false);
        
        bool isGlobal  = muon.isGlobalMuon();
        bool isTracker = muon.isTrackerMuon();
        bool isGlobalOrTracker = (isGlobal || isTracker);
        
        double norm_chi2     = 9999999.0;
        int minTrackerLayers = 0;
        int minValidMuonHits = 0;
        double _ip = 9999999.0;
        double _ipz = 9999999.0;
        int minPixelHits = 0;
        if ( muon.globalTrack().isNonnull() && muon.globalTrack().isAvailable() ){
            norm_chi2        = muon.normChi2();
            minTrackerLayers = static_cast<int> (muon.track()->hitPattern().trackerLayersWithMeasurement());
            minValidMuonHits = static_cast<int> (muon.globalTrack()->hitPattern().numberOfValidMuonHits());
            _ip = fabs(muon.dB());
            _ipz = fabs(muon.edB());
            minPixelHits = muon.innerTrack()->hitPattern().numberOfValidPixelHits();
        }
        
        int minMatchedStations = muon.numberOfMatches();
        
        double chIso = muon.userIsolation(pat::PfChargedHadronIso);
        double nhIso = muon.userIsolation(pat::PfNeutralHadronIso);
        double gIso  = muon.userIsolation(pat::PfGammaIso);
        double puIso = muon.userIsolation(pat::PfPUChargedHadronIso);
        double pt    = muon.pt() ;
        //double pfIso = (chIso + nhIso + gIso) / pt;
        double pfIso = (chIso + std::max(0.,nhIso + gIso - 0.5*puIso))/pt;
        
        if ( isGlobal  || ignoreCut("GlobalMuon")  )  passCut(ret, "GlobalMuon" );
        //else std::cout<<"failed GlobalMuon"<<std::endl;
        if ( isTracker || ignoreCut("TrackerMuon")  )  passCut(ret, "TrackerMuon" );
        //else std::cout<<"failed TrackerMuon"<<std::endl;
        if ( isGlobalOrTracker || ignoreCut("GlobalOrTrackerMuon")  )  passCut(ret, "GlobalOrTrackerMuon" );
        //else std::cout<<"failed GlobalOrTrackerMuon"<<std::endl;
        if ( norm_chi2          <  cut(indexChi2_,   double()) || ignoreCut(indexChi2_)    ) passCut(ret, indexChi2_   );
        //else std::cout<<"failed Chi2"<<std::endl;
        if ( minTrackerLayers   >  cut(indexMinTrackerLayers_,int()) || ignoreCut(indexMinTrackerLayers_)) passCut(ret, indexMinTrackerLayers_  );
        //else std::cout<<"failed minTrackerLayers"<<std::endl;
        if ( minValidMuonHits   >  cut(indexminValidMuHits_,int()) || ignoreCut(indexminValidMuHits_)) passCut(ret, indexminValidMuHits_  );
        //else std::cout<<"failed minValidMuonHits"<<std::endl;
        if ( _ip                <  cut(indexMaxIp_,double()) || ignoreCut(indexMaxIp_)) passCut(ret, indexMaxIp_  );
        //else std::cout<<"failed maxIp"<<std::endl;
        if ( _ipz               <  cut(indexMaxZImpact_,double()) || ignoreCut(indexMaxZImpact_)) passCut(ret, indexMaxZImpact_  );
        //else std::cout<<"failed maxZImpact"<<std::endl;
        if ( minPixelHits       >  cut(indexPixHits_,int())    || ignoreCut(indexPixHits_))  passCut(ret, indexPixHits_); 
        //else std::cout<<"failed minPixelHits"<<std::endl;
        if ( minMatchedStations >  cut(indexStations_,int())  || ignoreCut(indexStations_))  passCut(ret, indexStations_); 
        //else std::cout<<"failed minMatchedStations"<<std::endl;
        if ( pfIso              <  cut(indexmaxPfRelIso_, double())  || ignoreCut(indexmaxPfRelIso_)  ) passCut(ret, indexmaxPfRelIso_ ); 
        //else std::cout<<"failed pfIso"<<std::endl;
        
        setIgnored(ret);
        
        return (bool)ret;
    }
    
    
    
private: // member variables
    
    Version_t version_;
    
    index_type indexChi2_;
    index_type indexMinTrackerLayers_;
    index_type indexminValidMuHits_;
    index_type indexMaxIp_;
    index_type indexMaxZImpact_;
    index_type indexPixHits_;
    index_type indexStations_;
    index_type indexmaxPfRelIso_;
    
};

#endif
