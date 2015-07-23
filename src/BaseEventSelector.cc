#include <math.h>

#include "LJMet/Com/interface/BaseEventSelector.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

using namespace std;

BaseEventSelector::BaseEventSelector():
mName(""),
mLegend("")
{
}

void BaseEventSelector::BeginJob(std::map<std::string, edm::ParameterSet const > par)
{
    std::string _key = "event_selector";
    bool _missing_config = false;
    if ( par.find(_key)!=par.end() ){
        if (par[_key].exists("isMc")) mbPar["isMc"] = par[_key].getParameter<bool> ("isMc");
        else mbPar["isMc"] = false;
        
        if (par[_key].exists("btagOP")) msPar["btagOP"] = par[_key].getParameter<std::string> ("btagOP");
        else msPar["btagOP"] = "CSVM";
        
        if (par[_key].exists("JECup")) mbPar["JECup"] = par[_key].getParameter<bool> ("JECup");
        else mbPar["JECup"] = false;
        if (par[_key].exists("JECdown")) mbPar["JECdown"] = par[_key].getParameter<bool> ("JECdown");
        else mbPar["JECdown"] = false;
        if (par[_key].exists("JERup")) mbPar["JERup"] = par[_key].getParameter<bool> ("JERup");
        else mbPar["JERup"] = false;
        if (par[_key].exists("JERdown")) mbPar["JERdown"] = par[_key].getParameter<bool> ("JERdown");
        else mbPar["JERdown"] = false;
        if (par[_key].exists("JEC_txtfile")) msPar["JEC_txtfile"] = par[_key].getParameter<std::string> ("JEC_txtfile");
        else{
            msPar["JEC_txtfile"] = "";
            _missing_config = true;
        }
        if (par[_key].exists("BTagUncertUp")) mbPar["BTagUncertUp"] = par[_key].getParameter<bool> ("BTagUncertUp");
        else mbPar["BTagUncertUp"] = false;
        if (par[_key].exists("BTagUncertDown")) mbPar["BTagUncertDown"] = par[_key].getParameter<bool> ("BTagUncertDown");
        else mbPar["BTagUncertDown"] = false;
        
        if (par[_key].exists("MCL1JetPar")) msPar["MCL1JetPar"] = par[_key].getParameter<std::string> ("MCL1JetPar");
        else{
            msPar["MCL1JetPar"] = "../data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("MCL2JetPar")) msPar["MCL2JetPar"] = par[_key].getParameter<std::string> ("MCL2JetPar");
        else{
            msPar["MCL2JetPar"] = "../data/PHYS14_25_V2_L2Relative_AK4PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("MCL3JetPar")) msPar["MCL3JetPar"] = par[_key].getParameter<std::string> ("MCL3JetPar");
        else{
            msPar["MCL3JetPar"] = "../data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("MCL1JetParAK8")) msPar["MCL1JetParAK8"] = par[_key].getParameter<std::string> ("MCL1JetParAK8");
        else{
            msPar["MCL1JetParAK8"] = "../data/PHYS14_25_V2_L1FastJet_AK8PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("MCL2JetParAK8")) msPar["MCL2JetParAK8"] = par[_key].getParameter<std::string> ("MCL2JetParAK8");
        else{
            msPar["MCL2JetParAK8"] = "../data/PHYS14_25_V2_L2Relative_AK8PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("MCL3JetParAK8")) msPar["MCL3JetParAK8"] = par[_key].getParameter<std::string> ("MCL3JetParAK8");
        else{
            msPar["MCL3JetParAK8"] = "../data/PHYS14_25_V2_L3Absolute_AK8PFchs.txt";
            _missing_config = true;
        }
        
        if (par[_key].exists("DataL1JetPar")) msPar["DataL1JetPar"] = par[_key].getParameter<std::string> ("DataL1JetPar");
        else{
            msPar["DataL1JetPar"] = "../data/FT_53_V10_AN3_L1FastJet_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("DataL2JetPar")) msPar["DataL2JetPar"] = par[_key].getParameter<std::string> ("DataL2JetPar");
        else{
            msPar["DataL2JetPar"] = "../data/FT_53_V10_AN3_L2Relative_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("DataL3JetPar")) msPar["DataL3JetPar"] = par[_key].getParameter<std::string> ("DataL3JetPar");
        else{
            msPar["DataL3JetPar"] = "../data/FT_53_V10_AN3_L3Absolute_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("DataResJetPar")) msPar["DataResJetPar"] = par[_key].getParameter<std::string> ("DataResJetPar");
        else{
            msPar["DataResJetPar"] = "../data/FT_53_V10_AN3_L2L3Residual_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("DataL1JetParAK8")) msPar["DataL1JetParAK8"] = par[_key].getParameter<std::string> ("DataL1JetParAK8");
        else{
            msPar["DataL1JetParAK8"] = "../data/FT_53_V10_AN3_L1FastJet_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("DataL2JetParAK8")) msPar["DataL2JetParAK8"] = par[_key].getParameter<std::string> ("DataL2JetParAK8");
        else{
            msPar["DataL2JetParAK8"] = "../data/FT_53_V10_AN3_L2Relative_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("DataL3JetParAK8")) msPar["DataL3JetParAK8"] = par[_key].getParameter<std::string> ("DataL3JetParAK8");
        else{
            msPar["DataL3JetParAK8"] = "../data/FT_53_V10_AN3_L3Absolute_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("DataResJetParAK8")) msPar["DataResJetParAK8"] = par[_key].getParameter<std::string> ("DataResJetParAK8");
        else{
            msPar["DataResJetParAK8"] = "../data/FT_53_V10_AN3_L2L3Residual_AK5PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("doNewJEC")) mbPar["doNewJEC"] = par[_key].getParameter<bool> ("doNewJEC");
        else mbPar["doNewJEC"] = false;
        
        if (_missing_config) {
            std::cout << mLegend
            << "ONE OF THE FOLLOWING CONFIG OPTIONS MISSING!\n"
            << "MCL1JetPar, MCL2JetPar, MCL3JetPar, DataL1JetPar, DataL2JetPar,\n"
            << "DataL3JetPar, DataResJetPar" <<std::endl;
            std::cout << mLegend
            << "USING DEFAULT VALUES" << std::endl;
        }
    }
    
    msPar["btagger"] = mBtagCond.getAlgoName(msPar["btagOP"]);
    mdPar["btag_min_discr"] = mBtagCond.getDiscriminant(msPar["btagOP"]);
    
    bTagCut = mdPar["btag_min_discr"];
    std::cout << "b-tag check "<<msPar["btagOP"]<<" "<< msPar["btagger"]<<" "<<mdPar["btag_min_discr"]<<std::endl;
    
    if ( mbPar["isMc"] && ( mbPar["JECup"] || mbPar["JECdown"]))
        jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters(msPar["JEC_txtfile"].c_str(), "Total")));

    vector<JetCorrectorParameters> vPar;
    vector<JetCorrectorParameters> vParAK8;

    if ( mbPar["isMc"] && mbPar["doNewJEC"] ) {
        // Create the JetCorrectorParameter objects, the order does not matter.

        JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(msPar["MCL3JetPar"]);
        JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(msPar["MCL2JetPar"]);
    	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(msPar["MCL1JetPar"]);
        
	JetCorrectorParameters *L3JetParAK8  = new JetCorrectorParameters(msPar["MCL3JetParAK8"]);
        JetCorrectorParameters *L2JetParAK8  = new JetCorrectorParameters(msPar["MCL2JetParAK8"]);
    	JetCorrectorParameters *L1JetParAK8  = new JetCorrectorParameters(msPar["MCL1JetParAK8"]);
    	// Load the JetCorrectorParameter objects into a vector,
    	// IMPORTANT: THE ORDER MATTERS HERE !!!! 
    	vPar.push_back(*L1JetPar);
    	vPar.push_back(*L2JetPar);
    	vPar.push_back(*L3JetPar);

    	vParAK8.push_back(*L1JetParAK8);
    	vParAK8.push_back(*L2JetParAK8);
    	vParAK8.push_back(*L3JetParAK8);

    	std::cout << mLegend << "Applying new jet energy corrections" << std::endl;
    }
    else if ( !mbPar["isMc"] && mbPar["doNewJEC"] ) {
        // Create the JetCorrectorParameter objects, the order does not matter.

        JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(msPar["DataResJetPar"]); 
    	JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(msPar["DataL3JetPar"]);
    	JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(msPar["DataL2JetPar"]);
    	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(msPar["DataL1JetPar"]);

        JetCorrectorParameters *ResJetParAK8 = new JetCorrectorParameters(msPar["DataResJetParAK8"]); 
    	JetCorrectorParameters *L3JetParAK8  = new JetCorrectorParameters(msPar["DataL3JetParAK8"]);
    	JetCorrectorParameters *L2JetParAK8  = new JetCorrectorParameters(msPar["DataL2JetParAK8"]);
    	JetCorrectorParameters *L1JetParAK8  = new JetCorrectorParameters(msPar["DataL1JetParAK8"]);
    	// Load the JetCorrectorParameter objects into a vector,
    	// IMPORTANT: THE ORDER MATTERS HERE !!!! 
    	vPar.push_back(*L1JetPar);
   	vPar.push_back(*L2JetPar);
    	vPar.push_back(*L3JetPar);
    	vPar.push_back(*ResJetPar);

    	vParAK8.push_back(*L1JetParAK8);
   	vParAK8.push_back(*L2JetParAK8);
    	vParAK8.push_back(*L3JetParAK8);
    	vParAK8.push_back(*ResJetParAK8);

    	std::cout << mLegend << "Applying new jet energy corrections" << std::endl;
    }
    else{
    	std::cout << mLegend << "NOT applying new jet energy corrections - ARE YOU SURE?" << std::endl;
    }
    JetCorrector = new FactorizedJetCorrector(vPar);
    JetCorrectorAK8 = new FactorizedJetCorrector(vParAK8);
  
}

double BaseEventSelector::GetPerp(TVector3 & v1, TVector3 & v2)
{
    double perp;
    double _mag = v1.Cross(v2.Unit()).Mag();
    double _phi1 = v1.Phi();
    double _phi2 = v2.Phi();
    double _dphi = _phi1 - _phi2;
    if ( (_dphi > M_PI) || (_dphi > -M_PI && _dphi < 0.0) ) perp = _mag;
    else perp = -_mag;
    
    return perp;
}

void BaseEventSelector::Init( void )
{
    // init sanity check histograms
    mpEc->SetHistogram(mName, "jes_correction", 100, 0.8, 1.2);
    mpEc->SetHistogram(mName, "met_correction", 100, 0.0, 2.0);
    mpEc->SetHistogram(mName, "nBtagSfCorrections", 100, 0.0, 10.0);
}

TLorentzVector BaseEventSelector::correctJet(const pat::Jet & jet, edm::EventBase const & event, bool doAK8Corr)
{

  // JES and JES systematics
    pat::Jet correctedJet;
    if (mbPar["doNewJEC"])
        correctedJet = jet.correctedJet(0);                 //copy original jet
    else
        correctedJet = jet;                                 //copy default corrected jet

    double ptscale = 1.0;
    double unc = 1.0;
    double pt = correctedJet.pt();
    double correction = 1.0;

    edm::Handle<double> rhoHandle;
    edm::InputTag rhoSrc_("fixedGridRhoAll", "");
    event.getByLabel(rhoSrc_, rhoHandle);
    double rho = std::max(*(rhoHandle.product()), 0.0);

    if ( mbPar["isMc"] ){ 

    	if (mbPar["doNewJEC"]) {
      	    // We need to undo the default corrections and then apply the new ones

      	    double pt_raw = jet.correctedJet(0).pt();

	    if (doAK8Corr){
                JetCorrectorAK8->setJetEta(jet.eta());
          	JetCorrectorAK8->setJetPt(pt_raw);
                JetCorrectorAK8->setJetA(jet.jetArea());
          	JetCorrectorAK8->setRho(rho); 
    
                try{
    		    correction = JetCorrectorAK8->getCorrection();
                }
          	catch(...){
    		    std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
    		    std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
    		    std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
                }
	    }

	    else{
                JetCorrector->setJetEta(jet.eta());
          	JetCorrector->setJetPt(pt_raw);
                JetCorrector->setJetA(jet.jetArea());
          	JetCorrector->setRho(rho); 
    
                try{
    		    correction = JetCorrector->getCorrection();
                }
          	catch(...){
    		    std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
    		    std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
    		    std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
                }
	    }
      
            correctedJet.scaleEnergy(correction);
            pt = correctedJet.pt();

        }
        double factor = 0.0; // For Nominal Case
        double theAbsJetEta = fabs(jet.eta());
        
        if ( theAbsJetEta < 0.5 ) {
            factor = .052;
            if (mbPar["JERup"]) factor = 0.115;
            if (mbPar["JERdown"]) factor = -0.011;
        }
        else if ( theAbsJetEta < 1.1) {
            factor = 0.057;
            if (mbPar["JERup"]) factor = 0.114;
            if (mbPar["JERdown"]) factor = 0.0;
        }
        else if ( theAbsJetEta < 1.7) {
            factor = 0.096;
            if (mbPar["JERup"]) factor = 0.161;
            if (mbPar["JERdown"]) factor = 0.031;
        }
        else if ( theAbsJetEta < 2.3) {
            factor = 0.134;
            if (mbPar["JERup"]) factor = 0.228;
            if (mbPar["JERdown"]) factor = 0.040;
            
        }
        else if (theAbsJetEta < 5.0) {
            factor = 0.288;
            if (mbPar["JERup"]) factor = 0.488;
            if (mbPar["JERdown"]) factor = 0.088;
        }

        const reco::GenJet * genJet = jet.genJet();
    	if (genJet && genJet->pt()>15. && (fabs(genJet->pt()/pt-1)<0.5)){
      	    double gen_pt = genJet->pt();
      	    double reco_pt = pt;
            double deltapt = (reco_pt - gen_pt) * factor;
            ptscale = max(0.0, (reco_pt + deltapt) / reco_pt);
        }

        if ( mbPar["JECup"] || mbPar["JECdown"]) {
            jecUnc->setJetEta(jet.eta());
            jecUnc->setJetPt(pt*ptscale);

        if (mbPar["JECup"]) { 
	    try{
                unc = jecUnc->getUncertainty(true);
	    }
	    catch(...){ // catch all exceptions. Jet Uncertainty tool throws when binning out of range
	        std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
                std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
	        std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
	        unc = 0.0;
	    }
            unc = 1 + unc; 
        }
        else { 
	    try{
                unc = jecUnc->getUncertainty(false);
	    }
	    catch(...){
	        std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
	        std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
	        std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
	        unc = 0.0;
	    }
            unc = 1 - unc; 
        }

        if (pt*ptscale < 10.0 && mbPar["JECup"]) unc = 2.0;
        if (pt*ptscale < 10.0 && mbPar["JECdown"]) unc = 0.01;

        }
    }
    else if (!mbPar["isMc"]) {
      
        if (mbPar["doNewJEC"]) {
	
            // We need to undo the default corrections and then apply the new ones

            double pt_raw = jet.correctedJet(0).pt();
            JetCorrector->setJetEta(jet.eta());
            JetCorrector->setJetPt(pt_raw);
            JetCorrector->setJetA(jet.jetArea());
            JetCorrector->setRho(rho); 
	
	    try{
	        correction = JetCorrector->getCorrection();
	    }
	    catch(...){
	        std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
	        std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
	        std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
	    }
	  
            correctedJet.scaleEnergy(correction);
            pt = correctedJet.pt();

        }
    }

    TLorentzVector jetP4;
    jetP4.SetPtEtaPhiM(correctedJet.pt()*unc*ptscale, correctedJet.eta(),correctedJet.phi(), correctedJet.mass() );
    //std::cout<<"jet pt: "<<jetP4.Pt()<<" eta: "<<jetP4.Eta()<<" phi: "<<jetP4.Phi()<<" energy: "<<jetP4.E()<<std::endl;


    // sanity check - save correction of the first jet
    if (mNCorrJets==0){
        double _orig_pt = jet.pt();
        if (fabs(_orig_pt)<0.000000001){
            _orig_pt = 0.000000001;
        }
        SetHistValue("jes_correction", jetP4.Pt()/_orig_pt);
        ++mNCorrJets;
    }


    return jetP4;
}

bool BaseEventSelector::isJetTagged(const pat::Jet & jet, edm::EventBase const & event, bool applySF)
{
    bool _isTagged = false;
    
    if ( jet.bDiscriminator( msPar["btagger"] ) > bTagCut ) _isTagged = true;
    
    if (mbPar["isMc"] && applySF) {
        TLorentzVector lvjet = correctJet(jet, event);
        
        double _lightSf = mBtagCond.GetMistagScaleFactor(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        if ( mbPar["BTagUncertUp"] ) _lightSf += mBtagCond.GetMistagSFUncertUp(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        else if ( mbPar["BTagUncertDown"] ) _lightSf -= mBtagCond.GetMistagSFUncertDown(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        double _lightEff = mBtagCond.GetMistagRate(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        
        int _jetFlavor = abs(jet.partonFlavour());
        double _btagSf = mBtagCond.GetBtagScaleFactor(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        if ( mbPar["BTagUncertUp"] ) _btagSf += (mBtagCond.GetBtagSFUncertUp(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*(_jetFlavor==4?2:1));
        else if ( mbPar["BTagUncertDown"] ) _btagSf -= (mBtagCond.GetBtagSFUncertDown(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*(_jetFlavor==4?2:1));
        double _btagEff = mBtagCond.GetBtagEfficiency(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        
        mBtagSfUtil.SetSeed(abs(static_cast<int>(sin(jet.phi())*1e5)));
        
        // sanity check
        bool _orig_tag = _isTagged;
        
        mBtagSfUtil.modifyBTagsWithSF(_isTagged, _jetFlavor, _btagSf, _btagEff, _lightSf, _lightEff);
        
        // sanity check
        if (_isTagged != _orig_tag) ++mNBtagSfCorrJets;
        
    } // end of btag scale factor corrections
    return _isTagged;
}

TLorentzVector BaseEventSelector::correctMet(const pat::MET & met, edm::EventBase const & event)
{
    double correctedMET_px = met.px();
    double correctedMET_py = met.py();
    
    for (std::vector<edm::Ptr<pat::Jet> >::const_iterator ijet = mvAllJets.begin();
         ijet != mvAllJets.end(); ++ijet) {
        TLorentzVector lv = correctJet(**ijet, event);
        correctedMET_px += (**ijet).px() -lv.Px();
        correctedMET_py += (**ijet).py() -lv.Py();
    }
    
    correctedMET_p4.SetPxPyPzE(correctedMET_px, correctedMET_py, 0, sqrt(correctedMET_px*correctedMET_px+correctedMET_py*correctedMET_py));
    
    // sanity check histogram
    double _orig_met = met.pt();
    if (fabs(_orig_met) < 1.e-9) {
        _orig_met = 1.e-9;
    }
    SetHistValue("met_correction", correctedMET_p4.Pt()/_orig_met);
    return correctedMET_p4;
}
