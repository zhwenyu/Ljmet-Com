/*
 Contact:        Dylan Rankin (drankin@bu.edu)
 */
#ifndef LJMet_Com_interface_MVAElectronSelector_h
#define LJMet_Com_interface_MVAElectronSelector_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//Math
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "TMath.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBDT.h"

typedef math::XYZPoint Point;

class MVAElectronSelector : public Selector<pat::Electron>  {
    
public: // interface
    
    bool verbose_;
    
    enum Version_t { MEDIUM, TIGHT, NONE, N_VERSIONS};
    MVAElectronSelector() {}
    
    
    MVAElectronSelector( edm::ParameterSet const & parameters ){
        
        verbose_ = true;
        
        std::string versionStr = parameters.getParameter<std::string>("version");
        Version_t version = N_VERSIONS;
        
        if (versionStr == "MEDIUM"){
            version = MEDIUM;
            if(verbose_) std::cout << "MVAElectronSelector: You have choosen version = MEDIUM " <<  std::endl;
        }
        
        if (versionStr == "TIGHT"){
            version = TIGHT;
            if(verbose_) std::cout << "MVAElectronSelector: You have choosen version = TIGHT " <<  std::endl;
        }
        
        if ( versionStr == "NONE" ){
            version = NONE;
            if(verbose_){ std::cout << "MVAElectronSelector: If you want to use version NONE "
                << "then make sure to provide the selection cuts by yourself " << std::endl;}
        }

        initialize( version,
                   parameters.getParameter<Double_t>("cutValue_Bin"),
                   parameters.getParameter<Double_t>("cutValue_Bout"),
                   parameters.getParameter<Double_t>("cutValue_E")
                   );
        
        if ( parameters.exists("cutsToIgnore") )
            setIgnoredCuts( parameters.getParameter<std::vector<std::string> >("cutsToIgnore") );
       
        if ( parameters.exists("weightFiles") )
            weightFiles = parameters.getParameter<std::vector<std::string> >("weightFiles");
        if (weightFiles.size()!=3) {
            weightFiles.clear();
            weightFiles.push_back("../weights/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml");
            weightFiles.push_back("../weights/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml");
            weightFiles.push_back("../weights/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml");
        }
        // these are for 25ns, and are up-to-date as of Sep 24 2015
        // this needs to be checked periodically, as well as the list of variables for the MVA
        // look here: https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/RecoEgamma/ElectronIdentification/plugins/ElectronMVAEstimatorRun2Spring15NonTrig.cc
        // and here: https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Spring15_25ns_nonTrig_V1_cff.py
        retInternal_ = getBitTemplate();
        pvSrc_  = parameters.getParameter<edm::InputTag>("pvSrc");
        rhoSrc_ = parameters.getParameter<edm::InputTag>("rhoSrc");
        
    }
    
    void initialize(Version_t version, Double_t cutValue_Bin, Double_t cutValue_Bout, Double_t cutValue_E)
    {
        version_ = version;
        
        push_back("cutValue_Bin" );
        push_back("cutValue_Bout" );
        push_back("cutValue_E" );
        

        if (version_ == NONE){
            set("cutValue_Bin",  cutValue_Bin);
            set("cutValue_Bout",  cutValue_Bout);
            set("cutValue_E",  cutValue_E);
        }
        
        if (version_ == MEDIUM) {
            set("cutValue_Bin",  0.913286);
            set("cutValue_Bout",  0.805013);
            set("cutValue_E",  0.358969);
        }
        

        if (version_ == TIGHT) {
            set("cutValue_Bin",  0.967083);
            set("cutValue_Bout",  0.929117);
            set("cutValue_E",  0.726311);
        }
        
        indexCutValue_Bin_    = index_type(&bits_, "cutValue_Bin"   );
        indexCutValue_Bout_    = index_type(&bits_, "cutValue_Bout"   );
        indexCutValue_E_    = index_type(&bits_, "cutValue_E"   );
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
        const reco::BeamSpot &beamspot = *bsHandle.product();

        std::string mvaSpring15NonTrigWeightFile_V1;
        if (fabs(electron.superCluster()->eta())<=0.8 ) mvaSpring15NonTrigWeightFile_V1 = weightFiles[0];
        else if (fabs(electron.superCluster()->eta())<=1.479 ) mvaSpring15NonTrigWeightFile_V1 = weightFiles[1];
        else mvaSpring15NonTrigWeightFile_V1 = weightFiles[2];
       
        TMVA::Reader tmpTMVAReader( "!Color:Silent:!Error" );
      
        //
        // Configure all variables and spectators. Note: the order and names
        // must match what is found in the xml weights file!

        Float_t see, spp, OneMinusE1x5E5x5, R9, etawidth, phiwidth, HoE, PreShowerOverRaw, kfhits, kfchi2, gsfchi2, fbrem, convVtxFitProbability, EoP, eleEoPout, IoEmIoP, deta, dphi, detacalo, gsfhits, expectedMissingInnerHits, pt, isBarrel, isEndcap, SCeta, eClass, pfRelIso, expectedInnerHits, vtxconv, mcEventWeight, mcCBmatchingCategory;

        //
        // Pure ECAL -> shower shapes
        tmpTMVAReader.AddVariable("ele_oldsigmaietaieta", &see);
        tmpTMVAReader.AddVariable("ele_oldsigmaiphiiphi", &spp);
        tmpTMVAReader.AddVariable("ele_oldcircularity",   &OneMinusE1x5E5x5);
        tmpTMVAReader.AddVariable("ele_oldr9",            &R9);
        tmpTMVAReader.AddVariable("ele_scletawidth",      &etawidth);
        tmpTMVAReader.AddVariable("ele_sclphiwidth",      &phiwidth);
        tmpTMVAReader.AddVariable("ele_he",               &HoE);
        // Endcap only variables
        if( fabs(electron.superCluster()->eta())>1.479 )
          tmpTMVAReader.AddVariable("ele_psEoverEraw",    &PreShowerOverRaw);
        
        //Pure tracking variables
        tmpTMVAReader.AddVariable("ele_kfhits",           &kfhits);
        tmpTMVAReader.AddVariable("ele_kfchi2",           &kfchi2);
        tmpTMVAReader.AddVariable("ele_gsfchi2",        &gsfchi2);
      
        // Energy matching
        tmpTMVAReader.AddVariable("ele_fbrem",           &fbrem);
      
        tmpTMVAReader.AddVariable("ele_gsfhits",         &gsfhits);
        tmpTMVAReader.AddVariable("ele_expected_inner_hits",             &expectedMissingInnerHits);
        tmpTMVAReader.AddVariable("ele_conversionVertexFitProbability",  &convVtxFitProbability);
      
        tmpTMVAReader.AddVariable("ele_ep",              &EoP);
        tmpTMVAReader.AddVariable("ele_eelepout",        &eleEoPout);
        tmpTMVAReader.AddVariable("ele_IoEmIop",         &IoEmIoP);
        
        // Geometrical matchings
        tmpTMVAReader.AddVariable("ele_deltaetain",      &deta);
        tmpTMVAReader.AddVariable("ele_deltaphiin",      &dphi);
        tmpTMVAReader.AddVariable("ele_deltaetaseed",    &detacalo);

        // Spectator variables  
        tmpTMVAReader.AddSpectator("ele_pT",             &pt);
        tmpTMVAReader.AddSpectator("ele_isbarrel",       &isBarrel);
        tmpTMVAReader.AddSpectator("ele_isendcap",       &isEndcap);
        tmpTMVAReader.AddSpectator("scl_eta",            &SCeta);
      
        tmpTMVAReader.AddSpectator("ele_eClass",                 &eClass);
        tmpTMVAReader.AddSpectator("ele_pfRelIso",               &pfRelIso);
        tmpTMVAReader.AddSpectator("ele_expected_inner_hits",    &expectedInnerHits);
        tmpTMVAReader.AddSpectator("ele_vtxconv",                &vtxconv);
        tmpTMVAReader.AddSpectator("mc_event_weight",            &mcEventWeight);
        tmpTMVAReader.AddSpectator("mc_ele_CBmatching_category", &mcCBmatchingCategory);

        tmpTMVAReader.BookMVA( "Spring15_V1", mvaSpring15NonTrigWeightFile_V1 );

        // Pure ECAL -> shower shapes
        see            = electron.full5x5_sigmaIetaIeta();
        spp            = electron.full5x5_sigmaIphiIphi();
        OneMinusE1x5E5x5 = 1. - electron.full5x5_e1x5() / electron.full5x5_e5x5();
        R9             = electron.full5x5_r9();
        etawidth       = electron.superCluster()->etaWidth();
        phiwidth       = electron.superCluster()->phiWidth();
        HoE            = electron.hadronicOverEm();
        // Endcap only variables
        PreShowerOverRaw  = electron.superCluster()->preshowerEnergy() / electron.superCluster()->rawEnergy();
      
        // To get to CTF track information in pat::Electron, we have to have the pointer
        // to pat::Electron, it is not accessible from the pointer to reco::GsfElectron.
        // This behavior is reported and is expected to change in the future (post-7.4.5 some time).
        bool validKF= false; 
        reco::TrackRef myTrackRef = electron.closestCtfTrackRef();
        validKF = (myTrackRef.isAvailable() && (myTrackRef.isNonnull()) );  
      
        //Pure tracking variables
        kfhits         = (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;
        kfchi2          = (validKF) ? myTrackRef->normalizedChi2() : 0;
        gsfchi2         = electron.gsfTrack()->normalizedChi2();
      
        // Energy matching
        fbrem           = electron.fbrem();
      
        gsfhits         = electron.gsfTrack()->found();
        expectedMissingInnerHits = electron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
      
        reco::ConversionRef conv_ref = ConversionTools::matchedConversion(electron, conversions, beamspot.position());
        double vertexFitProbability = -1.; 
        if(!conv_ref.isNull()) {
          const reco::Vertex &vtx = conv_ref.get()->conversionVertex(); if (vtx.isValid()) {
            vertexFitProbability = TMath::Prob( vtx.chi2(), vtx.ndof());
          } 
        }
        convVtxFitProbability    = vertexFitProbability;
      
        EoP             = electron.eSuperClusterOverP();
        eleEoPout       = electron.eEleClusterOverPout();
        IoEmIoP         = (1.0/electron.ecalEnergy()) - (1.0 / electron.p());
      
        // Geometrical matchings
        deta            = electron.deltaEtaSuperClusterTrackAtVtx();
        dphi            = electron.deltaPhiSuperClusterTrackAtVtx();
        detacalo        = electron.deltaEtaSeedClusterTrackAtCalo();

        // Spectator variables  
        pt              = electron.pt();
        SCeta           = electron.superCluster()->eta();
        constexpr float ebeeSplit = 1.479;
        isBarrel        = ( fabs(SCeta) < ebeeSplit );
        isEndcap        = ( fabs(SCeta) >= ebeeSplit );
        // The spectator variables below were examined for training, but
        // are not necessary for evaluating the discriminator, so they are
        // given dummy values (the specator variables above are also unimportant).
        // They are introduced only to match the definition of the discriminator 
        // in the weights file.
        constexpr unsigned nines = 999;
        eClass               = nines;
        pfRelIso             = nines;
        expectedInnerHits    = nines;
        vtxconv              = nines;
        mcEventWeight        = nines;
        mcCBmatchingCategory = nines;

        // Constrain values

        if(fbrem < -1.) fbrem = -1.;
        deta = fabs(deta);
        if(deta > 0.06) deta = 0.06;
        dphi = fabs(dphi);
        if(dphi > 0.6) dphi = 0.6;
        if(EoP > 20.) EoP = 20.;
        if(eleEoPout > 20.) eleEoPout = 20.;
        detacalo = fabs(detacalo);
        if(detacalo > 0.2) detacalo = 0.2;
        if(OneMinusE1x5E5x5 < -1.) OneMinusE1x5E5x5 = -1;
        if(OneMinusE1x5E5x5 > 2.) OneMinusE1x5E5x5 = 2.; 
        if(R9 > 5) R9 = 5;
        if(gsfchi2 > 200.) gsfchi2 = 200;
        if(kfchi2 > 10.) kfchi2 = 10.;

        cutValue = tmpTMVAReader.EvaluateMVA( "Spring15_V1" );
	//std::cout<<"cutValue = "<<cutValue<<std::endl;

        return operator()(electron, ret);
    }
    
    //Spring 12 cuts
    bool operator()( const pat::Electron & electron, pat::strbitset & ret)
    {
        ret.set(false);
        Double_t scEta = electron.superCluster()->eta();
        

        // now apply the cuts
        if (fabs(scEta)<=0.8) {
            if (cutValue     >=  cut(indexCutValue_Bin_, double()) || ignoreCut(indexCutValue_Bin_) ) passCut(ret, indexCutValue_Bin_);
            passCut(ret, indexCutValue_Bout_);
            passCut(ret, indexCutValue_E_);
        }
        else if (fabs(scEta)<=1.479) {
            if (cutValue     >=  cut(indexCutValue_Bout_, double()) || ignoreCut(indexCutValue_Bout_) ) passCut(ret, indexCutValue_Bout_);
            passCut(ret, indexCutValue_Bin_);
            passCut(ret, indexCutValue_E_);
        }
        else {
            if (cutValue     >=  cut(indexCutValue_E_, double()) || ignoreCut(indexCutValue_E_) ) passCut(ret, indexCutValue_E_);
            passCut(ret, indexCutValue_Bin_);
            passCut(ret, indexCutValue_Bout_);
        }
        
        setIgnored(ret);

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
    double cutValue;
    std::vector<std::string> weightFiles;
    index_type indexCutValue_Bin_;
    index_type indexCutValue_Bout_;
    index_type indexCutValue_E_;
};

#endif
