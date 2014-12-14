/*
 Interface class for FWLite PAT analyzer-selectors
 Specific selectors must implement the () operator
 
 Author: Gena Kukartsev, 2010,2012
 */

#include <math.h>

#include "LJMet/Com/interface/BaseEventSelector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

using namespace std;

BaseEventSelector::BaseEventSelector():
mName(""),
mLegend("")
{
}

BaseEventSelector::~BaseEventSelector()
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
            msPar["MCL1JetPar"] = "../data/MCRUN2_72_V3A_L1FastJet_AK4PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("MCL2JetPar")) msPar["MCL2JetPar"] = par[_key].getParameter<std::string> ("MCL2JetPar");
        else{
            msPar["MCL2JetPar"] = "../data/MCRUN2_72_V3A_L2Relative_AK4PFchs.txt";
            _missing_config = true;
        }
        if (par[_key].exists("MCL3JetPar")) msPar["MCL3JetPar"] = par[_key].getParameter<std::string> ("MCL3JetPar");
        else{
            msPar["MCL3JetPar"] = "../data/MCRUN2_72_V3A_L3Absolute_AK4PFchs.txt";
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
}



void BaseEventSelector::EndJob()
{
}

void BaseEventSelector::AnalyzeEvent( edm::EventBase const & event, LjmetEventContent & ec )
{
}

std::string BaseEventSelector::GetName()
{
    return mName;
}

void BaseEventSelector::SetHistogram(std::string name, int nbins, double low, double high)
{
    // Declare a new histogram to be created for the module
    mpEc->SetHistogram(mName, name, nbins, low, high);
    return;
}

void BaseEventSelector::SetHistValue(std::string name,
                                     double value){
    mpEc->SetHistValue(mName, name, value);
    
    return;
}



// evaluates a signed perp components of v1 relative to v2
// the sign is defined by Phi
double BaseEventSelector::GetPerp(TVector3 & v1, TVector3 & v2){
    
    double perp;
    
    double _mag = v1.Cross(v2.Unit()).Mag();
    double _phi1 = v1.Phi();
    double _phi2 = v2.Phi();
    double _dphi = _phi1 - _phi2;
    if (_dphi>M_PI || (_dphi>-M_PI && _dphi<0)) perp = _mag;
    else perp = -_mag;
    
    return perp;
}



bool BaseEventSelector::AreMatched ( const reco::Candidate & c1,
                                    const reco::Candidate & c2,
                                    double DR,
                                    double DPtRel ) {
    unsigned int nPass=0;
    if ( (reco::deltaR(c1, c2) < DR) && (fabs(c2.pt() - c1.pt())/ c2.pt()<DPtRel)) {
        nPass++;
    }
    return (nPass>0);
}



std::vector<edm::Ptr<pat::Jet> > const & BaseEventSelector::GetAllJets() const { return mvAllJets; }
std::vector<edm::Ptr<pat::Jet> > const & BaseEventSelector::GetSelectedJets() const { return mvSelJets; }
std::vector<edm::Ptr<pat::Jet> > const & BaseEventSelector::GetLooseJets() const { return mvSelJets; }
std::vector<edm::Ptr<pat::Jet> > const & BaseEventSelector::GetSelectedBtagJets() const { return mvSelBtagJets; }
std::vector<std::pair<TLorentzVector,bool>> const & BaseEventSelector::GetCorrJetsWithBTags() const { return mvCorrJetsWithBTags; }
std::vector<edm::Ptr<pat::Muon> > const & BaseEventSelector::GetAllMuons() const { return mvAllMuons; }
std::vector<edm::Ptr<pat::Muon> > const & BaseEventSelector::GetSelectedMuons() const { return mvSelMuons; }
std::vector<edm::Ptr<pat::Muon> > const & BaseEventSelector::GetLooseMuons() const { return mvLooseMuons; }
std::vector<edm::Ptr<pat::Electron> > const & BaseEventSelector::GetAllElectrons() const { return mvAllElectrons; }
std::vector<edm::Ptr<pat::Electron> > const & BaseEventSelector::GetSelectedElectrons() const { return mvSelElectrons; }
std::vector<edm::Ptr<pat::Electron> > const & BaseEventSelector::GetLooseElectrons() const { return mvLooseElectrons; }
edm::Ptr<pat::MET> const & BaseEventSelector::GetMet() const { return mpMet; }
edm::Ptr<reco::PFMET> const & BaseEventSelector::GetType1CorrMet() const { return mpType1CorrMet; }
TLorentzVector const & BaseEventSelector::GetCorrectedMet() const { return correctedMET_p4; }
std::vector<unsigned > const & BaseEventSelector::GetSelectedTriggers() const { return mvSelTriggers; }
std::vector<edm::Ptr<reco::Vertex> > const & BaseEventSelector::GetSelectedPVs() const { return mvSelPVs; }

double const &
BaseEventSelector::GetTestValue() const{
    return mTestValue;
}



void BaseEventSelector::SetTestValue(double & test){
    mTestValue = test;
    return;
}



void BaseEventSelector::SetCorrJetsWithBTags(std::vector<std::pair<TLorentzVector,bool>> & jets){
    mvCorrJetsWithBTags = jets;
    return;
}



void BaseEventSelector::SetCorrectedMet(TLorentzVector & met){
    correctedMET_p4 = met;
    return;
}



void BaseEventSelector::SetMc(bool isMc){
    mbIsMc = isMc;
    return;
}



bool BaseEventSelector::IsMc(){
    return mbIsMc;
}

bool BaseEventSelector::isJetTagged(const pat::Jet & jet, edm::EventBase const & event, bool applySF)
{
    bool _isTagged = false;
    
    
    if ( jet.bDiscriminator( msPar["btagger"] ) > bTagCut ) _isTagged = true;
    // DEBUG: list all available taggers with values
    // cout << "_isTagged "<<_isTagged<<endl;
    
    // const std::vector<std::pair<std::string, float> > & vpd = jet.getPairDiscri();
    // for ( std::vector<std::pair<std::string, float> >::const_iterator i = vpd.begin();
    // i != vpd.end(); ++i){
    // std::cout << mLegend << i->first << ", " << i->second << std::endl;
    // }
    
    if (mbPar["isMc"] && applySF){
        
        TLorentzVector lvjet = correctJet(jet, event);
        
        double _lightSf = mBtagCond.GetMistagScaleFactor(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        if ( mbPar["BTagUncertUp"] ) _lightSf += mBtagCond.GetMistagSFUncertUp(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        else if ( mbPar["BTagUncertDown"] )_lightSf -= mBtagCond.GetMistagSFUncertDown(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        double _lightEff = mBtagCond.GetMistagRate(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        
        int _jetFlavor = abs(jet.partonFlavour());
        double _btagSf = mBtagCond.GetBtagScaleFactor(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        if ( mbPar["BTagUncertUp"] ) _btagSf += (mBtagCond.GetBtagSFUncertUp(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*(_jetFlavor==4?2:1));
        else if ( mbPar["BTagUncertDown"] )_btagSf -= (mBtagCond.GetBtagSFUncertDown(lvjet.Et(), lvjet.Eta(), msPar["btagOP"])*(_jetFlavor==4?2:1));
        double _btagEff = mBtagCond.GetBtagEfficiency(lvjet.Et(), lvjet.Eta(), msPar["btagOP"]);
        
        mBtagSfUtil.SetSeed(abs(static_cast<int>(sin(jet.phi())*100000)));
        
        // sanity check
        bool _orig_tag = _isTagged;
        
        mBtagSfUtil.modifyBTagsWithSF(_isTagged, _jetFlavor, _btagSf, _btagEff,
                                      _lightSf, _lightEff);
        
        // sanity check
        if (_isTagged != _orig_tag) ++mNBtagSfCorrJets;
        
    } // end of btag scale factor corrections
    return _isTagged;
    
}

TLorentzVector BaseEventSelector::correctJet(const pat::Jet & jet, edm::EventBase const & event)
{
    // JES and JES systematics
    pat::Jet correctedJet;
    correctedJet = jet; //copy original jet
    
    double ptscale = 1.0;
    double unc = 1.0;
    double pt = correctedJet.pt();
    double correction = 1.0;
    
    edm::Handle<double> rhoHandle;
    edm::InputTag rhoSrc_("fixedGridRhoAll", "");
    event.getByLabel(rhoSrc_, rhoHandle);
    double rho = std::max(*(rhoHandle.product()), 0.0);
    
    if ( mbPar["isMc"] ) {
        double factor = 0.0; // For Nominal Case
        double theAbsJetEta = abs(jet.eta());
        
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
        if (genJet && genJet->pt()>15. && (abs(genJet->pt()/pt-1)<0.5)) {
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
                try {
                    unc = jecUnc->getUncertainty(false);
                }
                catch(...) {
                    std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
                    std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
                    std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
                    unc = 0.0;
                }
                unc = 1 - unc;
            }
            
            if (pt*ptscale < 10.0) {
                if (mbPar["JECup"]) unc = 2.0;
                if (mbPar["JECdown"]) unc = 0.01;
            }
        }
    }
    
    TLorentzVector jetP4;
    jetP4.SetPtEtaPhiM(correctedJet.pt()*unc*ptscale, correctedJet.eta(),correctedJet.phi(), correctedJet.mass() );
    
    // sanity check - save correction of the first jet
    if (mNCorrJets == 0) {
        double _orig_pt = jet.pt();
        if (fabs(_orig_pt)<1.e-9) {
            _orig_pt = 1.e-9;
        }
        SetHistValue("jes_correction", jetP4.Pt()/_orig_pt);
        ++mNCorrJets;
    }
    return jetP4;
}




// Start of public:

void BaseEventSelector::Init( void )
{
    // init sanity check histograms
    mpEc->SetHistogram(mName, "jes_correction", 100, 0.8, 1.2);
    mpEc->SetHistogram(mName, "met_correction", 100, 0.0, 2.0);
    mpEc->SetHistogram(mName, "nBtagSfCorrections", 100, 0.0, 10.0);
    
    return;
}

void BaseEventSelector::SetEventContent(LjmetEventContent * pEc)
{
    mpEc = pEc;
    return;
}


// Start of protected:

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


// Start of private:

void BaseEventSelector::init()
{
    // private init method to be called by LjmetFactory when registering the selector
    mLegend = "[" + mName + "]: ";
    std::cout << mLegend << "registering " << mName << std::endl;
    return;
}

void BaseEventSelector::setName(std::string name)
{
    mName = name;
}

void BaseEventSelector::BeginEvent(edm::EventBase const & event, LjmetEventContent & ec )
{
    // Do what any event selector must do before event gets checked
    mNCorrJets = 0;
    mNBtagSfCorrJets = 0;
    return;
}

void BaseEventSelector::EndEvent(edm::EventBase const & event, LjmetEventContent & ec )
{
    // Do what any event selector must do after event processing is done
    // (but before event content gets saved to file)
    SetHistValue("nBtagSfCorrections", mNBtagSfCorrJets);
    return;
}
