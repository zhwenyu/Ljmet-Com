/*
 Calculator for the SLiTT analysis
 
 Author: Gena Kukartsev, 2012
 */



#include <iostream>
#include <string>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "LJMet/Com/interface/TopElectronSelector.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "LJMet/Com/interface/MiniIsolation.h"
#include "HepPDT/ParticleID.hh"


using std::cout;
using std::endl;

class LjmetFactory;



class DileptonCalc : public BaseCalc {
    
public:
    
    DileptonCalc();
    virtual ~DileptonCalc(){}
    
    virtual int BeginJob();
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}
    
private:
    
    bool                      isMc;
    std::string               dataType;
    edm::InputTag             rhoSrc_it;
    edm::InputTag             pvCollection_it;
    edm::InputTag             genParticles_it;
    edm::InputTag             packedPFCandsLabel_;
    edm::InputTag             triggerObjects_;
    edm::InputTag             triggerBits_;
    std::vector<unsigned int> keepPDGID;
    std::vector<unsigned int> keepMomPDGID;
    bool keepFullMChistory;
    bool doTriggerStudy_; 
    double rhoIso;

    boost::shared_ptr<TopElectronSelector>     electronSelL_, electronSelM_, electronSelT_;
    std::vector<reco::Vertex> goodPVs;
    int findMatch(const reco::GenParticleCollection & genParticles, int idToMatch, double eta, double phi);
    double mdeltaR(double eta1, double phi1, double eta2, double phi2);
    void fillMotherInfo(const reco::Candidate *mother, int i, std::vector <int> & momid, std::vector <int> & momstatus, std::vector<double> & mompt, std::vector<double> & mometa, std::vector<double> & momphi, std::vector<double> & momenergy);
};

static int reg = LjmetFactory::GetInstance()->Register(new DileptonCalc(), "DileptonCalc");

DileptonCalc::DileptonCalc()
{
}

int DileptonCalc::BeginJob()
{
    if (mPset.exists("dataType"))     dataType = mPset.getParameter<std::string>("dataType");
    else                              dataType = "None";
    
    if (mPset.exists("rhoSrc"))       rhoSrc_it = mPset.getParameter<edm::InputTag>("rhoSrc");
    else                              rhoSrc_it = edm::InputTag("fixedGridRhoFastjetAll", "", "");
    
    if (mPset.exists("pvCollection")) pvCollection_it = mPset.getParameter<edm::InputTag>("pvCollection");
    else                              pvCollection_it = edm::InputTag("offlineSlimmedPrimaryVertices");
    
    if (mPset.exists("isMc"))         isMc = mPset.getParameter<bool>("isMc");
    else                              isMc = false;
    
    if (mPset.exists("genParticles")) genParticles_it = mPset.getParameter<edm::InputTag>("genParticles");
    else                              genParticles_it = edm::InputTag("prunedGenParticles");

    if (mPset.exists("packedPFCands"))	packedPFCandsLabel_ = mPset.getParameter<edm::InputTag>("packedPFCands");
    else                              	packedPFCandsLabel_ = edm::InputTag("packedPFCandidates");

    if (mPset.exists("keepPDGID"))    keepPDGID = mPset.getParameter<std::vector<unsigned int> >("keepPDGID");
    else                              keepPDGID.clear();
    
    if (mPset.exists("keepMomPDGID")) keepMomPDGID = mPset.getParameter<std::vector<unsigned int> >("keepMomPDGID");
    else                              keepMomPDGID.clear();

    if(mPset.exists("TriggerObjects")) triggerObjects_ = mPset.getParameter<edm::InputTag>("TriggerObjects");
    else                               triggerObjects_ = edm::InputTag("selectedPatTrigger");

    if(mPset.exists("TriggerBits")) triggerBits_ = mPset.getParameter<edm::InputTag>("TriggerBits");
    else                            triggerBits_ = edm::InputTag("TriggerResults","","HLT");
    
    if(mPset.exists("doTriggerStudy")) doTriggerStudy_ = mPset.getParameter<bool>("doTriggerStudy");
    else                             doTriggerStudy_ = true;

    //   if (mPset.exists("keepFullMChistory")) keepFullMChistory = mPset.getParameter<bool>("keepFullMChistory");
    //else                                   keepFullMChistory = true;
    keepFullMChistory = true;
    cout << "keepFullMChistory "     <<    keepFullMChistory << endl;
    
    if ( mPset.exists("cutbasedIDSelectorLoose")){
        electronSelL_ = boost::shared_ptr<TopElectronSelector>(
                                                               new TopElectronSelector(mPset.getParameter<edm::ParameterSet>("cutbasedIDSelectorLoose")) );
    }
    else {
        std::cout << "DileptonCalc: Loose electron selector not configured, exiting"
        << std::endl;
        std::exit(-1);
    }
    if ( mPset.exists("cutbasedIDSelectorMedium")){
        electronSelM_ = boost::shared_ptr<TopElectronSelector>(
                                                               new TopElectronSelector(mPset.getParameter<edm::ParameterSet>("cutbasedIDSelectorMedium")) );
    }
    else {
        std::cout << "DileptonCalc: Medium electron selector not configured, exiting"
        << std::endl;
        std::exit(-1);
    }
    if ( mPset.exists("cutbasedIDSelectorTight")){
        electronSelT_ = boost::shared_ptr<TopElectronSelector>(
                                                               new TopElectronSelector(mPset.getParameter<edm::ParameterSet>("cutbasedIDSelectorTight")) );
    }
    else {
        std::cout << "DileptonCalc: Tight electron selector not configured, exiting"
        << std::endl;
        std::exit(-1);
    }
    
   
    return 0;
}

int DileptonCalc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector)
{
    //
    // compute event variables here
    //
    
    //
    // _____ Get objects from the selector _____________________
    //
    std::vector<edm::Ptr<pat::Muon> >     const & vSelMuons        = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons    = selector->GetSelectedElectrons();
    std::vector<edm::Ptr<pat::Jet> >      const & vSelJets         = selector->GetSelectedJets();
    std::vector<pat::Jet>                 const & vSelCleanedJets  = selector->GetSelectedCleanedJets();
    edm::Ptr<pat::MET>                    const & pMet             = selector->GetMet();
    std::vector<unsigned int>             const & vSelTriggers     = selector->GetSelectedTriggers();
    
    //
    // _____ Primary dataset (from python cfg) _____________________
    //
    //
    int dataEE = 0;
    int dataEM = 0;
    int dataMM = 0;
    
    if      (dataType == "EE" or dataType == "ElEl") dataEE = 1;
    else if (dataType == "EM" or dataType == "ElMu") dataEM = 1;
    else if (dataType == "MM" or dataType == "MuMu") dataMM = 1;
    else if (dataType == "All" or dataType == "ALL") {
        dataEE = 1; dataEM = 1; dataMM = 1;
    }
    
    SetValue("dataEE", dataEE);
    SetValue("dataEM", dataEM);
    SetValue("dataMM", dataMM);
    
    //
    // ____ Trigger ____________________________
    //
    int passEE = 0;
    int passEM = 0;
    int passMM = 0;
    
    if (vSelTriggers.size() == 3) {
        passEE = (int)vSelTriggers.at(0);
        passEM = (int)vSelTriggers.at(1);
        passMM = (int)vSelTriggers.at(2);
    }
    
    SetValue("trigEE", passEE);
    SetValue("trigEM", passEM);
    SetValue("trigMM", passMM);


    bool HLT_DoubleEle33=false;
    bool HLT_DoubleEle33_MW=false;
    bool HLT_DoubleEle8_Mass8_HT300 = false;
    bool HLT_Ele17_Ele12_DZ=false;
    bool HLT_Ele27WP85=false;
    bool HLT_Mu27TkMu8=false;
    bool HLT_Mu30TkMu11=false;
    bool HLT_Mu40TkMu11=false;
    bool HLT_Mu17_Iso_Mu8_Iso=false;
    bool HLT_Mu17_Iso_TkMu8_Iso=false;
    bool HLT_DoubleMu8_Mass8_HT300 = false;
    bool HLT_Mu40=false;
    bool HLT_IsoTkMu24=false;
    bool HLT_DoubleMu33NoFiltersNoVtx=false;
    bool HLT_Mu17Ele12=false;
    bool HLT_Mu8Ele17=false;
    bool HLT_Mu23Ele12=false;
    bool HLT_Mu8Ele23=false;
    bool HLT_Mu8_Ele8_Mass8_HT300 = false;
    bool HLT_Mu30Ele30=false;
    bool HLT_PFHT900=false;
    bool HLT_AK8PFJet360TrimMass30=false;

    if(doTriggerStudy_){

      edm::Handle<edm::TriggerResults> triggerBits;
      event.getByLabel(triggerBits_,triggerBits);
      const edm::TriggerNames &names = event.triggerNames(*triggerBits);
      
      
      for (unsigned int i=0; i!=triggerBits->size(); ++i) {
	std::string Path = names.triggerName(i);
	
	const unsigned int triggerIndex(i);
	if(triggerBits->accept(triggerIndex)){
	  //electron paths
	  if(Path=="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1" || Path=="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2" || Path=="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3") HLT_DoubleEle33=true;
	  if(Path=="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v1"|| Path=="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v2" || Path=="HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3") HLT_DoubleEle33_MW=true;
	  if(Path=="HLT_Ele27_eta2p1_WP85_Gsf_v1" || Path=="HLT_Ele27_eta2p1_WP85_Gsf_v2" || Path=="HLT_Ele27_eta2p1_WP85_Gsf_v3") HLT_Ele27WP85=true;
	  if(Path=="HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1" || Path=="HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2" || Path=="HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3")HLT_Ele17_Ele12_DZ=true;
	  if(Path=="HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1" || Path=="HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v2" || Path=="HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v3") HLT_DoubleEle8_Mass8_HT300 = true;
	  //Muon paths
	  if(Path=="HLT_Mu27_TkMu8_v1" || Path=="HLT_Mu27_TkMu8_v2" || Path=="HLT_Mu27_TkMu8_v3") HLT_Mu27TkMu8=true;
	  if(Path=="HLT_Mu30_TkMu11_v1" || Path=="HLT_Mu30_TkMu11_v2" || Path=="HLT_Mu30_TkMu11_v3") HLT_Mu30TkMu11=true;
	  if(Path=="HLT_Mu40_TkMu11_v1" || Path=="HLT_Mu40_TkMu11_v2" || Path=="HLT_Mu40_TkMu11_v3") HLT_Mu40TkMu11=true;
	  if(Path=="HLT_DoubleMu8_CaloIdM_TrackIdM_Mass8_PFHT300_v1" || Path=="HLT_DoubleMu8_CaloIdM_TrackIdM_Mass8_PFHT300_v2" || Path=="HLT_DoubleMu8_CaloIdM_TrackIdM_Mass8_PFHT300_v3") HLT_DoubleMu8_Mass8_HT300 = true;
	  if(Path=="HLT_M40_v1" || Path=="HLT_M40_v2" || Path=="HLT_M40_v3") HLT_Mu40=true;
	  if(Path=="HLT_IsoTkMu24_IterTrk02_v1" || Path=="HLT_IsoTkMu24_IterTrk02_v2" || Path=="HLT_IsoTkMu24_IterTrk02_v3") HLT_IsoTkMu24=true;
	  if(Path=="HLT_DoubleMu33NoFiltersNoVtx_v1" || Path=="HLT_DoubleMu33NoFiltersNoVtx_v2" || Path=="HLT_DoubleMu33NoFiltersNoVtx_v3") HLT_DoubleMu33NoFiltersNoVtx=true;
	  if(Path=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1" || Path=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2" || Path=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3") HLT_Mu17_Iso_Mu8_Iso=true;
	  if(Path=="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1" || Path=="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2" || Path=="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3") HLT_Mu17_Iso_TkMu8_Iso=true;
	  //cross paths
	  if(Path=="HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1" || Path=="HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2" || Path=="HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3") HLT_Mu17Ele12=true;
	  if(Path=="HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1" || Path=="HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2" || Path=="HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3") HLT_Mu8Ele17=true;
	  if(Path=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1" || Path=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2" || Path=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3") HLT_Mu23Ele12=true;
	  if(Path=="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1" || Path=="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2" || Path=="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3") HLT_Mu8Ele23=true;
	  if(Path=="HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v1" || Path=="HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v2" || Path=="HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v3") HLT_Mu30Ele30=true;
	  if(Path=="HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1" || Path=="HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v2" || Path=="HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v3") HLT_Mu8_Ele8_Mass8_HT300 = true;
	  //HT/Jet
	  if(Path=="HLT_PFHT900_v1" || Path=="HLT_PFHT900_v2" || Path=="HLT_PFHT900_v3") HLT_PFHT900=true;
	  if(Path=="HLT_AK8PFJet360TrimMod_Mass30_v1" || Path=="HLT_AK8PFJet360TrimMod_Mass30_v2" || Path=="HLT_AK8PFJet360TrimMod_Mass30_v3") HLT_AK8PFJet360TrimMass30=true;
	}
      }
    }


    SetValue("HLT_DoubleEle33",HLT_DoubleEle33);
    SetValue("HLT_DoubleEle33_MW",HLT_DoubleEle33_MW);
    SetValue("HLT_DoubleEle8_Mass8_HT300",HLT_DoubleEle8_Mass8_HT300);
    SetValue("HLT_Ele17_Ele12_DZ",HLT_Ele17_Ele12_DZ);
    SetValue("HLT_Ele27WP85",HLT_Ele27WP85);
    SetValue("HLT_Mu27TkMu8",HLT_Mu27TkMu8);
    SetValue("HLT_Mu30TkMu11",HLT_Mu30TkMu11);
    SetValue("HLT_Mu40TkMu11",HLT_Mu40TkMu11);
    SetValue("HLT_Mu17_Iso_Mu8_Iso",HLT_Mu17_Iso_Mu8_Iso);
    SetValue("HLT_Mu17_Iso_TkMu8_Iso",HLT_Mu17_Iso_TkMu8_Iso);
    SetValue("HLT_DoubleMu8_Mass8_HT300",HLT_DoubleMu8_Mass8_HT300);
    SetValue("HLT_Mu40",HLT_Mu40);
    SetValue("HLT_IsoTkMu24",HLT_IsoTkMu24);
    SetValue("HLT_DoubleMu33NoFiltersNoVtx",HLT_DoubleMu33NoFiltersNoVtx);
    SetValue("HLT_Mu17Ele12",HLT_Mu17Ele12);
    SetValue("HLT_Mu8Ele17",HLT_Mu8Ele17);
    SetValue("HLT_Mu23Ele12",HLT_Mu23Ele12);
    SetValue("HLT_Mu8Ele23",HLT_Mu8Ele23);
    SetValue("HLT_Mu30Ele30",HLT_Mu30Ele30);
    SetValue("HLT_Mu8_Ele8_Mass8_HT300",HLT_Mu8_Ele8_Mass8_HT300);
    SetValue("HLT_PFHT900",HLT_PFHT900);
    SetValue("HLT_AK8PFJet360TrimMass30",HLT_AK8PFJet360TrimMass30);


    //
    //_____ Event kinematics __________________
    //
    
    //Primary vertices
    edm::Handle<std::vector<reco::Vertex> > pvHandle;
    event.getByLabel(pvCollection_it, pvHandle);
    goodPVs = *(pvHandle.product());
    
    SetValue("nPV", (int)goodPVs.size());
    
    //
    //_____ Electrons _________________________
    //
    
    //Four std::vector
    std::vector <double> elPt;
    std::vector <double> elEta;
    std::vector <double> elPhi;
    std::vector <double> elEnergy;
    
    //Quality criteria
    std::vector <double> elRelIso;
    std::vector <double> elMiniIsoEA;
    std::vector <double> elMiniIsoDB;
    std::vector <double> elDxy;
    std::vector <int>    elNotConversion;
    std::vector <int>    elChargeConsistent;
    std::vector <int>    elIsEBEE;
    std::vector <int>    elQuality;
    std::vector <int>    elCharge;
    std::vector <int>    elGsfCharge;
    std::vector <int>    elCtfCharge;
    std::vector <int>    elScPixCharge;
    
    //ID requirement
    std::vector <double> elDeta;
    std::vector <double> elDphi;
    std::vector <double> elSihih;
    std::vector <double> elHoE;
    std::vector <double> elD0;
    std::vector <double> elSIP3D;
    std::vector <double> elDZ;
    std::vector <double> elOoemoop;
    std::vector <int>    elMHits;
    std::vector <int>    elVtxFitConv;
    std::vector<double>  elMVA;

    //added CMSDAS variables
    std::vector <double> diElMass;
    std::vector <int> elCharge1;
    std::vector <int> elCharge2;

    //Extra info about isolation
    std::vector <double> elChIso;
    std::vector <double> elNhIso;
    std::vector <double> elPhIso;
    std::vector <double> elAEff;
    std::vector <double> elRhoIso;
    std::vector <double> elPUIso;
    
    //mother-information
    //Generator level information -- MC matching
    std::vector<double> elGen_Reco_dr;
    std::vector<int> elPdgId;
    std::vector<int> elStatus;
    std::vector<int> elMatched;
    std::vector<int> elNumberOfMothers;
    std::vector<double> elMother_pt;
    std::vector<double> elMother_eta;
    std::vector<double> elMother_phi;
    std::vector<double> elMother_energy;
    std::vector<int> elMother_id;
    std::vector<int> elMother_status;
    //Matched gen electron information:
    std::vector<double> elMatchedPt;
    std::vector<double> elMatchedEta;
    std::vector<double> elMatchedPhi;
    std::vector<double> elMatchedEnergy;
    std::vector<std::string> TriggerElectronFilters;
    std::vector<double> TriggerElectronPts;
    std::vector<double> TriggerElectronEtas;
    std::vector<double> TriggerElectronPhis;
    std::vector<double> TriggerElectronEnergies;
    
    edm::Handle<double> rhoHandle;
    event.getByLabel(rhoSrc_it, rhoHandle);
    rhoIso = std::max(*(rhoHandle.product()), 0.0);

    //rho isolation from susy recommendation
    edm::Handle<double> rhoJetsNC;
    event.getByLabel(edm::InputTag("fixedGridRhoFastjetCentralNeutral","") , rhoJetsNC);
    double myRhoJetsNC = *rhoJetsNC;
    double _rhoNC = myRhoJetsNC;

    //pfcandidates for mini iso
    
    edm::Handle<pat::PackedCandidateCollection> packedPFCands;
    event.getByLabel(packedPFCandsLabel_, packedPFCands);
	
    
    pat::strbitset retElectron  = electronSelL_->getBitTemplate();
    bool retElectronT,retElectronM,retElectronL;
    

    
    //
    //_____Electrons______
    //
    
    //keep track of which electron we are looking at
    int ElIndex = 0;
    for (std::vector<edm::Ptr<pat::Electron> >::const_iterator iel = vSelElectrons.begin(); iel != vSelElectrons.end(); iel++){
      //Protect against electrons without tracks (should never happen, but just in case)
      if ((*iel)->gsfTrack().isNonnull() and (*iel)->gsfTrack().isAvailable()){
	
	//Four std::vector
	elPt     . push_back((*iel)->ecalDrivenMomentum().pt()); //Must check: why ecalDrivenMomentum?
	elEta    . push_back((*iel)->ecalDrivenMomentum().eta());
	elPhi    . push_back((*iel)->ecalDrivenMomentum().phi());
	elEnergy . push_back((*iel)->ecalDrivenMomentum().energy());
        
	//if there are two electrons calculate invariant mass from two highest pt objects
	if(vSelElectrons.size()==2){
	  for (std::vector<edm::Ptr<pat::Electron> >::const_iterator iiel = vSelElectrons.begin(); iiel != vSelElectrons.end(); iiel++){
	    //float mass = pow( (*iel)->energy() * (*iiel)->energy() - ( (*iel)->pt() * (*iiel)->pt() - (*iel)->pz() * (*iiel)->pz()), 0.5);
	    if(iiel!=iel){
	      TLorentzVector diElFourVec( (*iel)->px() + (*iiel)->px(),(*iel)->py() + (*iiel)->py(),(*iel)->pz() + (*iiel)->pz(), (*iel)->energy() + (*iiel)->energy());
	      diElMass.push_back(diElFourVec.M());
	      elCharge1.push_back( (*iel)->charge());
	      elCharge2.push_back( (*iiel)->charge());
	    }
	  }
	}
	else{ 
	  diElMass.push_back(-1);
	  elCharge1.push_back(-99999);
	  elCharge2.push_back(-99998);
	}
	//Isolation
	//double AEff  = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03,
	//							       (*iel)->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);

	//implement effective area: up-to-date as of PHYS14
	double AEff;
	if( fabs((*iel)->ecalDrivenMomentum().eta())<1.0) AEff=0.1752;
	else if(fabs((*iel)->ecalDrivenMomentum().eta())<1.479) AEff=0.1862;
	else if(fabs((*iel)->ecalDrivenMomentum().eta())<2.0) AEff=0.1411;
	else if(fabs((*iel)->ecalDrivenMomentum().eta())<2.2) AEff=0.1534;
	else if(fabs((*iel)->ecalDrivenMomentum().eta())<2.3) AEff=0.1903;
	else if(fabs((*iel)->ecalDrivenMomentum().eta())<2.4) AEff=0.2243;
	else if(fabs((*iel)->ecalDrivenMomentum().eta())<2.5) AEff=0.2687;
	
	reco::GsfElectron::PflowIsolationVariables pfIso = (*iel)->pfIsolationVariables();
	double chIso = pfIso.sumChargedHadronPt;
	double nhIso = pfIso.sumNeutralHadronEt;
	double phIso = pfIso.sumPhotonEt;
	double PUIso = pfIso.sumPUPt;
	double relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) ) / (*iel)->pt();
	
	elChIso  . push_back(chIso);
	elNhIso  . push_back(nhIso);
	elPhIso  . push_back(phIso);
	elPUIso  . push_back(PUIso);
	elAEff   . push_back(AEff);
	elRhoIso . push_back(rhoIso);
        
	elRelIso . push_back(relIso);
        
	//Conversion rejection
	int nLostHits = (*iel)->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
	double dist   = (*iel)->convDist();
	double dcot   = (*iel)->convDcot();
	int notConv   = nLostHits == 0 and (fabs(dist) > 0.02 or fabs(dcot) > 0.02);
	//get three different charges
	elGsfCharge.push_back( (*iel)->gsfTrack()->charge());
	if( (*iel)->closestCtfTrackRef().isNonnull()) elCtfCharge.push_back((*iel)->closestCtfTrackRef()->charge());
	else elCtfCharge.push_back(-999);
	elScPixCharge.push_back((*iel)->scPixCharge());
	elCharge.push_back((*iel)->charge());
	elNotConversion . push_back(notConv);
        
	retElectronL = (*electronSelL_)(**iel, event, retElectron);
	retElectronM = (*electronSelM_)(**iel, event, retElectron);
	retElectronT = (*electronSelT_)(**iel, event, retElectron);
        
	elQuality.push_back((retElectronT<<2) + (retElectronM<<1) + retElectronL);
        
	//IP: for some reason this is with respect to the first vertex in the collection
	if(goodPVs.size() > 0){
	  elDxy.push_back((*iel)->gsfTrack()->dxy(goodPVs.at(0).position()));
	  elDZ.push_back((*iel)->gsfTrack()->dz(goodPVs.at(0).position()));
	} else {
	  elDxy.push_back(-999);
	  elDZ.push_back(-999);
	}
	elChargeConsistent.push_back((*iel)->isGsfCtfScPixChargeConsistent());
	elIsEBEE.push_back(((*iel)->isEBEEGap()<<2) + ((*iel)->isEE()<<1) + (*iel)->isEB());
	elDeta.push_back((*iel)->deltaEtaSuperClusterTrackAtVtx());
	elDphi.push_back((*iel)->deltaPhiSuperClusterTrackAtVtx());
	elSihih.push_back((*iel)->full5x5_sigmaIetaIeta());
	elHoE.push_back((*iel)->hcalOverEcal());
	elD0.push_back((*iel)->dB());
	pat::Electron::IpType elIP3d = (pat::Electron::IpType) 1;
	float sip3d = (*iel)->dB(elIP3d) / (*iel)->edB(elIP3d);
	elSIP3D.push_back(sip3d);
	elOoemoop.push_back(1.0/(*iel)->ecalEnergy() - (*iel)->eSuperClusterOverP()/(*iel)->ecalEnergy());
	elMHits.push_back((*iel)->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
	elVtxFitConv.push_back((*iel)->passConversionVeto());
        //add mva
	elMVA.push_back( selector->mvaValue( *(iel->get()), event) );

	//add miniIso
	elMiniIsoEA.push_back(getPFMiniIsolation_EffectiveArea(packedPFCands, dynamic_cast<const reco::Candidate *>(iel->get()),0.05, 0.2, 10., false, false,myRhoJetsNC));
	elMiniIsoDB.push_back(getPFMiniIsolation_DeltaBeta(packedPFCands, dynamic_cast<const reco::Candidate *>(iel->get()), 0.05, 0.2, 10., false));


	//Trigger Matching - store 4-std::vector and filter information for all trigger objects deltaR matched to electrons
	/*if(doTriggerStudy_){
	  
	  //read in trigger objects
	  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	  event.getByLabel(triggerObjects_,triggerObjects);
	  
	  edm::Handle<edm::TriggerResults> triggerBits;
	  event.getByLabel(triggerBits_,triggerBits);
	  const edm::TriggerNames &names = event.triggerNames(*triggerBits);
	  
	  //loop over them for deltaR matching
	  TLorentzVector trigObj;
	  pat::TriggerObjectStandAlone matchedObj;
	  std::vector<std::string> paths;
	  float closestDR = 10000.;
	  for( pat::TriggerObjectStandAlone obj : *triggerObjects){
	    obj.unpackPathNames(names);
	    float dR = mdeltaR( (*iel)->eta(), (*iel)->phi(), obj.eta(),obj.phi() );
	    if(dR < closestDR){
	      closestDR = dR;
	      trigObj.SetPtEtaPhiE((*iel)->pt(),(*iel)->eta(),(*iel)->phi(),(*iel)->energy());
	      matchedObj=obj;
	    }
	  }
	  if(closestDR<0.5){
	    TriggerElectronPts.push_back(trigObj.Pt());
	    TriggerElectronEtas.push_back(trigObj.Eta());
	    TriggerElectronPhis.push_back(trigObj.Phi());
	    TriggerElectronEnergies.push_back(trigObj.Energy());
	    //std::cout<<"found matched trigger object!"<<std::endl;
	    //now store information about filters
	    for (unsigned h = 0; h < matchedObj.filterLabels().size(); ++h){
	      std::string filter = matchedObj.filterLabels()[h];
	      std::string Index = Form("MatchedIndex%i_",ElIndex);
	      std::string hltFilter_wIndex = Index+filter;
	      //std::cout<<hltFilter_wIndex<<std::endl;
	      TriggerElectronFilters.push_back(hltFilter_wIndex);
	    }	    
	  }
	  else{
	    TriggerElectronPts.push_back(-9999);
	    TriggerElectronEtas.push_back(-9999);
	    TriggerElectronPhis.push_back(-9999);
	    TriggerElectronEnergies.push_back(-9999);
	  }
	  }
	*/

	if(isMc && keepFullMChistory){
	  //cout << "start\n";
	  edm::Handle<reco::GenParticleCollection> genParticles;
	  event.getByLabel(genParticles_it, genParticles);
	  int matchId = findMatch(*genParticles, 11, (*iel)->eta(), (*iel)->phi());
	  double closestDR = 10000.;
	  //cout << "matchId "<<matchId <<endl;
	  if (matchId>=0) {
	    const reco::GenParticle & p = (*genParticles).at(matchId);
	    closestDR = mdeltaR( (*iel)->eta(), (*iel)->phi(), p.eta(), p.phi());
	    //cout << "closestDR "<<closestDR <<endl;
	    if(closestDR < 0.3){
	      elGen_Reco_dr.push_back(closestDR);
	      elPdgId.push_back(p.pdgId());
	      elStatus.push_back(p.status());
	      elMatched.push_back(1);
	      elMatchedPt.push_back( p.pt());
	      elMatchedEta.push_back(p.eta());
	      elMatchedPhi.push_back(p.phi());
	      elMatchedEnergy.push_back(p.energy());
	      int oldSize = elMother_id.size();
	      fillMotherInfo(p.mother(), 0, elMother_id, elMother_status, elMother_pt, elMother_eta, elMother_phi, elMother_energy);
	      elNumberOfMothers.push_back(elMother_id.size()-oldSize);
	    }
	  }
	  if(closestDR >= 0.3){
	    elNumberOfMothers.push_back(-1);
	    elGen_Reco_dr.push_back(-1.0);
	    elPdgId.push_back(-1);
	    elStatus.push_back(-1);
	    elMatched.push_back(0);
	    elMatchedPt.push_back(-1000.0);
	    elMatchedEta.push_back(-1000.0);
	    elMatchedPhi.push_back(-1000.0);
	    elMatchedEnergy.push_back(-1000.0);
            
	  }
          
          
	}//closing the isMC checking criteria
      }
      //increment index
      ElIndex+=1;
    }
    
    //trigger info
    SetValue("TrigElPt",TriggerElectronPts);
    SetValue("TrigElEta",TriggerElectronEtas);
    SetValue("TrigElPhi",TriggerElectronPhis);
    SetValue("TrigElEnergy",TriggerElectronEnergies);
    SetValue("TrigElFilters",TriggerElectronFilters);

    //Four std::vector
    SetValue("elPt"     , elPt);
    SetValue("elEta"    , elEta);
    SetValue("elPhi"    , elPhi);
    SetValue("elEnergy" , elEnergy);

    //cmsdas variables
    SetValue("diElMass", diElMass);
    SetValue("elCharge1", elCharge1);
    SetValue("elCharge2", elCharge2);
    
    SetValue("elGsfCharge",elGsfCharge);
    SetValue("elCtfCharge",elCtfCharge);
    SetValue("elScPixCharge",elScPixCharge);
    SetValue("elCharge", elCharge);
    //Quality requirements
    SetValue("elRelIso" , elRelIso); //Isolation
    SetValue("elDxy"    , elDxy);    //Dxy
    SetValue("elNotConversion" , elNotConversion);  //Conversion rejection
    SetValue("elChargeConsistent", elChargeConsistent);
    SetValue("elIsEBEE", elIsEBEE);
    SetValue("elQuality", elQuality);
    SetValue("elMiniIsoEA",elMiniIsoEA);
    SetValue("elMiniIsoDB",elMiniIsoDB);
    //this value not specific to electrons but we set it here
    SetValue("rhoNC",_rhoNC);

    //ID cuts
    SetValue("elDeta", elDeta);
    SetValue("elDphi", elDphi);
    SetValue("elSihih", elSihih);
    SetValue("elHoE", elHoE);
    SetValue("elD0", elD0);
    SetValue("elSIP3D",elSIP3D);
    SetValue("elDZ", elDZ);
    SetValue("elOoemoop", elOoemoop);
    SetValue("elMHits", elMHits);
    SetValue("elVtxFitConv", elVtxFitConv);
    SetValue("elMVA",elMVA);

    //Extra info about isolation
    SetValue("elChIso" , elChIso);
    SetValue("elNhIso" , elNhIso);
    SetValue("elPhIso" , elPhIso);
    SetValue("elAEff"  , elAEff);
    SetValue("elRhoIso", elRhoIso);
    SetValue("elPUIso", elPUIso);
    
    //MC matching -- mother information
    SetValue("elNumberOfMothers", elNumberOfMothers);
    SetValue("elGen_Reco_dr", elGen_Reco_dr);
    SetValue("elPdgId", elPdgId);
    SetValue("elStatus", elStatus);
    SetValue("elMatched",elMatched);
    SetValue("elMother_pt", elMother_pt);
    SetValue("elMother_eta", elMother_eta);
    SetValue("elMother_phi", elMother_phi);
    SetValue("elMother_energy", elMother_energy);
    SetValue("elMother_status", elMother_status);
    SetValue("elMother_id", elMother_id);
    //Matched gen muon information:
    SetValue("elMatchedPt", elMatchedPt);
    SetValue("elMatchedEta", elMatchedEta);
    SetValue("elMatchedPhi", elMatchedPhi);
    SetValue("elMatchedEnergy", elMatchedEnergy);
    
    
    //
    //_____ Muons _____________________________
    //
    
    std::vector <int> muCharge;
    std::vector <bool> muGlobal;
    std::vector <bool> muTracker;
    std::vector <bool> muPF;
    //Four std::vector
    std::vector <double> muPt;
    std::vector <double> muEta;
    std::vector <double> muPhi;
    std::vector <double> muEnergy;
    
    //Quality criteria
    std::vector <double> muChi2;
    std::vector <double> muDxy;
    std::vector <double> muDz;
    std::vector <double> muSIP3D;
    std::vector <double> muRelIso;
    std::vector <double> muMiniIsoEA;
    std::vector <double> muMiniIsoDB;
    
    std::vector <int> muNValMuHits;
    std::vector <int> muNMatchedStations;
    std::vector <int> muNValPixelHits;
    std::vector <int> muNTrackerLayers;
    
    //Extra info about isolation
    std::vector <double> muChIso;
    std::vector <double> muNhIso;
    std::vector <double> muGIso;
    std::vector <double> muPuIso;

    //ID info
    std::vector <int> muIsTight;
    std::vector<int> muIsLoose;
    
    //Generator level information -- MC matching
    std::vector<double> muGen_Reco_dr;
    std::vector<int> muPdgId;
    std::vector<int> muStatus;
    std::vector<int> muMatched;
    std::vector<int> muNumberOfMothers;
    std::vector<double> muMother_pt;
    std::vector<double> muMother_eta;
    std::vector<double> muMother_phi;
    std::vector<double> muMother_energy;
    std::vector<int> muMother_id;
    std::vector<int> muMother_status;
    //Matched gen muon information:
    std::vector<double> muMatchedPt;
    std::vector<double> muMatchedEta;
    std::vector<double> muMatchedPhi;
    std::vector<double> muMatchedEnergy;

    std::vector<std::string> TriggerMuonFilters;
    std::vector<double> TriggerMuonPts;
    std::vector<double> TriggerMuonEtas;
    std::vector<double> TriggerMuonPhis;
    std::vector<double> TriggerMuonEnergies;
    

    //make index for muons
    int MuIndex = 0;
    for (std::vector<edm::Ptr<pat::Muon> >::const_iterator imu = vSelMuons.begin(); imu != vSelMuons.end(); imu++){
        //Protect against muons without tracks (should never happen, but just in case)
        if ((*imu)->globalTrack().isNonnull() and (*imu)->globalTrack().isAvailable() and
            (*imu)->innerTrack().isNonnull()  and (*imu)->innerTrack().isAvailable()){
            
            
            //charge
            muCharge.push_back((*imu)->charge());
            
            //Four std::vector
            muPt     . push_back((*imu)->pt());
            muEta    . push_back((*imu)->eta());
            muPhi    . push_back((*imu)->phi());
            muEnergy . push_back((*imu)->energy());
            

	    muIsTight.push_back((*imu)->isTightMuon(goodPVs.at(0)));
	    muIsLoose.push_back((*imu)->isLooseMuon());

            muGlobal.push_back((*imu)->isGlobalMuon());
	    muTracker.push_back((*imu)->isTrackerMuon());
	    muPF.push_back((*imu)->isPFMuon());
            //Chi2
            muChi2 . push_back((*imu)->globalTrack()->normalizedChi2());
            
            //Isolation
	    /*            double chIso  = (*imu)->userIsolation(pat::PfChargedHadronIso);
            double nhIso  = (*imu)->userIsolation(pat::PfNeutralHadronIso);
            double gIso   = (*imu)->userIsolation(pat::PfGammaIso);
            double puIso  = (*imu)->userIsolation(pat::PfPUChargedHadronIso);*/

	    //new definition of iso based on muon pog page
	    const reco::MuonPFIsolation pfIsolationR04 = (*imu)->pfIsolationR04();
	    double chIso  = pfIsolationR04.sumChargedHadronPt;
            double nhIso  = pfIsolationR04.sumNeutralHadronEt;
            double gIso   = pfIsolationR04.sumPhotonEt;
            double puIso  = pfIsolationR04.sumPUPt;
            double relIso = (chIso + std::max(0.,nhIso + gIso - 0.5*puIso)) / (*imu)->pt();
            muRelIso . push_back(relIso);

            muChIso . push_back(chIso);
            muNhIso . push_back(nhIso);
            muGIso  . push_back(gIso);
            muPuIso . push_back(puIso);

            //get miniIso
	    muMiniIsoEA.push_back(getPFMiniIsolation_EffectiveArea(packedPFCands, dynamic_cast<const reco::Candidate *>(imu->get()),0.05, 0.2, 10., false, false,myRhoJetsNC));
	    muMiniIsoDB.push_back(getPFMiniIsolation_DeltaBeta(packedPFCands, dynamic_cast<const reco::Candidate *>(imu->get()), 0.05, 0.2, 10., false));            
	    float musip3d;
	    pat::Muon::IpType muIP3d = (pat::Muon::IpType) 1;
            //IP: for some reason this is with respect to the first vertex in the collection
            if (goodPVs.size() > 0){
	      muDxy . push_back((*imu)->dB());
	      muDz  . push_back((*imu)->muonBestTrack()->dz(goodPVs.at(0).position()));
	      //muon sip3d;
	      musip3d = (*imu)->dB(muIP3d) / (*imu)->edB(muIP3d);

            } 
	    else {
	      muDxy . push_back(-999);
	      muDz  . push_back(-999);
	      musip3d=-999;
            }
	    muSIP3D.push_back(musip3d);
            //Numbers of hits
            muNValMuHits       . push_back((*imu)->globalTrack()->hitPattern().numberOfValidMuonHits());
            muNMatchedStations . push_back((*imu)->numberOfMatchedStations());
            muNValPixelHits    . push_back((*imu)->innerTrack()->hitPattern().numberOfValidPixelHits());
            muNTrackerLayers   . push_back((*imu)->innerTrack()->hitPattern().trackerLayersWithMeasurement());

	   
	    //Trigger Matching - store 4-std::vector and filter information for all trigger objects deltaR matched to electrons
	    /*if(doTriggerStudy_){
	  
	      //read in trigger objects
	      edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	      event.getByLabel(triggerObjects_,triggerObjects);
	      
	      edm::Handle<edm::TriggerResults> triggerBits;
	      event.getByLabel(triggerBits_,triggerBits);
	      const edm::TriggerNames &names = event.triggerNames(*triggerBits);
	      
	      //loop over them for deltaR matching
	      TLorentzVector trigObj;
	      pat::TriggerObjectStandAlone matchedObj;
	      std::vector<std::string> paths;
	      float closestDR = 10000.;
	      for( pat::TriggerObjectStandAlone obj : *triggerObjects){
		obj.unpackPathNames(names);
		float dR = mdeltaR( (*imu)->eta(), (*imu)->phi(), obj.eta(),obj.phi() );
		if(dR < closestDR){
		  closestDR = dR;
		  trigObj.SetPtEtaPhiE((*imu)->pt(),(*imu)->eta(),(*imu)->phi(),(*imu)->energy());
		  matchedObj=obj;
		}
	      }
	      if(closestDR<0.5){
		TriggerMuonPts.push_back(trigObj.Pt());
		TriggerMuonEtas.push_back(trigObj.Eta());
		TriggerMuonPhis.push_back(trigObj.Phi());
		TriggerMuonEnergies.push_back(trigObj.Energy());
		//std::cout<<"found muon matched trigger object!"<<std::endl;
		//now store information about filters
		for (unsigned int h = 0; h < matchedObj.filterLabels().size(); h++){
		  std::string filter = matchedObj.filterLabels()[h];
		  std::string Index = Form("MatchedIndex%i_",MuIndex);
		  std::string hltFilter_wIndex = Index+filter;
		  //std::cout<<hltFilter_wIndex<<std::endl;
		  TriggerMuonFilters.push_back(hltFilter_wIndex);
		}	    
	      }
	      else{
		TriggerMuonPts.push_back(-9999);
		TriggerMuonEtas.push_back(-9999);
		TriggerMuonPhis.push_back(-9999);
		TriggerMuonEnergies.push_back(-9999);
	      }
	      
	    }*/

            
            if(isMc && keepFullMChistory){
                edm::Handle<reco::GenParticleCollection> genParticles;
                event.getByLabel(genParticles_it, genParticles);
                int matchId = findMatch(*genParticles, 13, (*imu)->eta(), (*imu)->phi());
                double closestDR = 10000.;
                if (matchId>=0) {
                    const reco::GenParticle & p = (*genParticles).at(matchId);
                    closestDR = mdeltaR( (*imu)->eta(), (*imu)->phi(), p.eta(), p.phi());
                    if(closestDR < 0.3){
                        muGen_Reco_dr.push_back(closestDR);
                        muPdgId.push_back(p.pdgId());
                        muStatus.push_back(p.status());
                        muMatched.push_back(1);
                        muMatchedPt.push_back( p.pt());
                        muMatchedEta.push_back(p.eta());
                        muMatchedPhi.push_back(p.phi());
                        muMatchedEnergy.push_back(p.energy());
                        int oldSize = muMother_id.size();
                        fillMotherInfo(p.mother(), 0, muMother_id, muMother_status, muMother_pt, muMother_eta, muMother_phi, muMother_energy);
                        muNumberOfMothers.push_back(muMother_id.size()-oldSize);
                    }
                }
                if(closestDR >= 0.3){
                    muNumberOfMothers.push_back(-1);
                    muGen_Reco_dr.push_back(-1.0);
                    muPdgId.push_back(-1);
                    muStatus.push_back(-1);
                    muMatched.push_back(0);
                    muMatchedPt.push_back(-1000.0);
                    muMatchedEta.push_back(-1000.0);
                    muMatchedPhi.push_back(-1000.0);
                    muMatchedEnergy.push_back(-1000.0);
                    
                }
                
                
            }
        }
	//increment index
	MuIndex+=1;
    }

    //trigger info
    SetValue("TrigMuPt",TriggerMuonPts);
    SetValue("TrigMuEta",TriggerMuonEtas);
    SetValue("TrigMuPhi",TriggerMuonPhis);
    SetValue("TrigMuEnergy",TriggerMuonEnergies);
    SetValue("TrigMuFilters",TriggerMuonFilters);
    
    
    SetValue("muCharge", muCharge);
    SetValue("muGlobal", muGlobal);
    SetValue("muTracker",muTracker);
    SetValue("muPF",muPF);

    //Four std::vector
    SetValue("muPt"     , muPt);
    SetValue("muEta"    , muEta);
    SetValue("muPhi"    , muPhi);
    SetValue("muEnergy" , muEnergy);
  
    //muon ids
    SetValue("muIsTight", muIsTight);
    SetValue("muIsLoose",muIsLoose);
  
    //Quality criteria
    SetValue("muChi2"   , muChi2);
    SetValue("muDxy"    , muDxy);
    SetValue("muDz"     , muDz);
    SetValue("muSIP3D", muSIP3D);
    SetValue("muRelIso" , muRelIso);
    SetValue("muMiniIsoEA", muMiniIsoEA);
    SetValue("muMiniIsoDB", muMiniIsoDB);

    SetValue("muNValMuHits"       , muNValMuHits);
    SetValue("muNMatchedStations" , muNMatchedStations);
    SetValue("muNValPixelHits"    , muNValPixelHits);
    SetValue("muNTrackerLayers"   , muNTrackerLayers);
    
    //Extra info about isolation
    SetValue("muChIso", muChIso);
    SetValue("muNhIso", muNhIso);
    SetValue("muGIso" , muGIso);
    SetValue("muPuIso", muPuIso);
    
    //MC matching -- mother information
    SetValue("muGen_Reco_dr", muGen_Reco_dr);
    SetValue("muPdgId", muPdgId);
    SetValue("muStatus", muStatus);
    SetValue("muMatched",muMatched);
    SetValue("muMother_pt", muMother_pt);
    SetValue("muMother_eta", muMother_eta);
    SetValue("muMother_phi", muMother_phi);
    SetValue("muMother_energy", muMother_energy);
    SetValue("muMother_status", muMother_status);
    SetValue("muMother_id", muMother_id);
    SetValue("muNumberOfMothers", muNumberOfMothers);
    //Matched gen muon information:
    SetValue("muMatchedPt", muMatchedPt);
    SetValue("muMatchedEta", muMatchedEta);
    SetValue("muMatchedPhi", muMatchedPhi);
    SetValue("muMatchedEnergy", muMatchedEnergy);


    //
    //_____GenJets_____________________________
    //

    if(isMc){
      edm::InputTag genJetColl = edm::InputTag("slimmedGenJets");
      edm::Handle<std::vector<reco::GenJet> > genJets;
      
      event.getByLabel(genJetColl,genJets);
    
      std::vector<double> genJetPt;
      std::vector<double> genJetEta;
      std::vector<double> genJetPhi;
      std::vector<double> genJetEnergy;
      std::vector<double> genJetMass;
      
      for(std::vector<reco::GenJet>::const_iterator ijet = genJets->begin(); ijet != genJets->end(); ijet++){
	
	genJetPt.push_back(ijet->pt());
	genJetEta.push_back(ijet->eta());
	genJetPhi.push_back(ijet->phi());
	genJetEnergy.push_back(ijet->energy());
	genJetMass.push_back(ijet->mass());
	//cout<<"setting gen jet info"<<endl;
	
      }
      
      SetValue("genJetPt", genJetPt);
      SetValue("genJetEta", genJetEta);
      SetValue("genJetPhi", genJetPhi);
      SetValue("genJetEnergy", genJetEnergy);
      SetValue("genJetMass", genJetMass);
    }

    //
    //_____ Jets ______________________________
    //
    
    //Get Top-like jets
    edm::InputTag topJetColl = edm::InputTag("slimmedJetsAK8");
    edm::Handle<std::vector<pat::Jet> > topJets;
    event.getByLabel(topJetColl, topJets);
    
    //Four std::vector -- COMMENT OUT BECAUSE FOR NOW WE CAN'T ACCESS THESE
    std::vector <double> CATopJetPt;
    std::vector <double> CATopJetEta;
    std::vector <double> CATopJetPhi;
    std::vector <double> CATopJetEnergy;
    
    std::vector <double> CATopJetCSV;
    //   std::vector <double> CATopJetRCN;
    
    //Identity
    std::vector <int> CATopJetIndex;
    std::vector <int> CATopJetnDaughters;
    
    //Top-like properties
    std::vector <double> CATopJetTopMass;
    std::vector <double> CATopJetMinPairMass;
    
    //Daughter four std::vector and index
    std::vector <double> CATopDaughterPt;
    std::vector <double> CATopDaughterEta;
    std::vector <double> CATopDaughterPhi;
    std::vector <double> CATopDaughterEnergy;
    
    std::vector <int> CATopDaughterMotherIndex;
    
    for (std::vector<pat::Jet>::const_iterator ijet = topJets->begin(); ijet != topJets->end(); ijet++) {
        
        int index = (int)(ijet-topJets->begin());
        
	CATopJetPt     . push_back(ijet->pt());
        CATopJetEta    . push_back(ijet->eta());
        CATopJetPhi    . push_back(ijet->phi());
        CATopJetEnergy . push_back(ijet->energy());
        
        CATopJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
        //     CATopJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));
	 
        CATopJetIndex      . push_back(index);
        CATopJetnDaughters . push_back((int)ijet->numberOfDaughters());
        
	//	cout<<"tag infos"<<(ijet->tagInfoLabels()).at(0)<<endl;

        if(ijet->tagInfoLabels().size()>0){
	  reco::CATopJetTagInfo const * jetInfo = dynamic_cast<reco::CATopJetTagInfo const *> (ijet->tagInfo("caTop"));
	  CATopJetTopMass     . push_back(jetInfo->properties().topMass);
	  CATopJetMinPairMass . push_back(jetInfo->properties().minMass);
	}
	else{
	  CATopJetTopMass     . push_back(-999);
	  CATopJetMinPairMass . push_back(-999);
	}

        for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++) {
	  //cout<<"in top daughter loop"<<endl;
	  CATopDaughterPt     . push_back(ijet->daughter(ui)->pt());
	  CATopDaughterEta    . push_back(ijet->daughter(ui)->eta());
	  CATopDaughterPhi    . push_back(ijet->daughter(ui)->phi());
	  CATopDaughterEnergy . push_back(ijet->daughter(ui)->energy());
	  
	  CATopDaughterMotherIndex . push_back(index);
        }
	//cout<<"finished top daughter loop"<<endl;
    }
    
    //Four std::vector
    SetValue("CATopJetPt"    , CATopJetPt);
    SetValue("CATopJetEta"   , CATopJetEta);
    SetValue("CATopJetPhi"   , CATopJetPhi);
    SetValue("CATopJetEnergy", CATopJetEnergy);
    
    SetValue("CATopJetCSV"   , CATopJetCSV);
    //   SetValue("CATopJetRCN"    , CATopJetRCN);
    
    //Identity
    SetValue("CATopJetIndex"      , CATopJetIndex);
    SetValue("CATopJetnDaughters" , CATopJetnDaughters);
    
    //Properties
    SetValue("CATopJetTopMass"     , CATopJetTopMass);
    SetValue("CATopJetMinPairMass" , CATopJetMinPairMass);
    //cout<<"setting values for top daughters"<<endl;
    //Daughter four std::vector and index
    SetValue("CATopDaughterPt"     , CATopDaughterPt);
    SetValue("CATopDaughterEta"    , CATopDaughterEta);
    SetValue("CATopDaughterPhi"    , CATopDaughterPhi);
    SetValue("CATopDaughterEnergy" , CATopDaughterEnergy);
    
    SetValue("CATopDaughterMotherIndex"      , CATopDaughterMotherIndex);
    
    //Get CA8 jets for W's
    edm::InputTag CAWJetColl = edm::InputTag("slimmedJetsAK8");
    edm::Handle<std::vector<pat::Jet> > CAWJets;
    event.getByLabel(CAWJetColl, CAWJets);
    
    //Four std::vector
    std::vector <double> CAWJetPt;
    std::vector <double> CAWJetEta;
    std::vector <double> CAWJetPhi;
    std::vector <double> CAWJetEnergy;
    
    std::vector <double> CAWJetCSV;
    //   std::vector <double> CAWJetRCN;
    
    //Identity
    std::vector <int> CAWJetIndex;
    std::vector <int> CAWJetnDaughters;
    
    //Mass
    std::vector <double> CAWJetTrimmedMass;
    std::vector <double> CAWJetPrunedMass;
    std::vector <double> CAWJetFilteredMass;
    
    //nsubjettiness
    std::vector<double> CAWJetTau1;
    std::vector<double> CAWJetTau2;
    std::vector<double> CAWJetTau3;


    //Daughter four std::vector and index -- THESE ARE CURRENTLY IDENTICAL TO THOSE FOR TOP DAUGHTERS BECAUSE THE JET IS ALWAYS THE SLIMMED AK8JET
    std::vector <double> CAWDaughterPt;
    std::vector <double> CAWDaughterEta;
    std::vector <double> CAWDaughterPhi;
    std::vector <double> CAWDaughterEnergy;
    
    std::vector <int> CAWDaughterMotherIndex;
    
    for (std::vector<pat::Jet>::const_iterator ijet = CAWJets->begin(); ijet != CAWJets->end(); ijet++){
        
        int index = (int)(ijet-CAWJets->begin());
        
        //Four std::vector
	CAWJetPt     . push_back(ijet->pt());
        CAWJetEta    . push_back(ijet->eta());
        CAWJetPhi    . push_back(ijet->phi());
        CAWJetEnergy . push_back(ijet->energy());
        
        CAWJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
        //     CAWJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));
        
        //Identity
        CAWJetIndex      . push_back(index);
        CAWJetnDaughters . push_back((int)ijet->numberOfDaughters());
	//cout<<"about to set w masses"<<endl;
        //Mass
        CAWJetTrimmedMass . push_back(ijet->userFloat("ak8PFJetsCHSTrimmedLinks"));
        CAWJetPrunedMass . push_back(ijet->userFloat("ak8PFJetsCHSPrunedLinks"));
        CAWJetFilteredMass . push_back(ijet->userFloat("ak8PFJetsCHSFilteredLinks"));
	//cout<<"set w masses"<<endl;
	//nsubjettiness
	CAWJetTau1.push_back( ijet->userFloat("NjettinessAK8:tau1"));
	CAWJetTau2.push_back( ijet->userFloat("NjettinessAK8:tau2"));
	CAWJetTau3.push_back( ijet->userFloat("NjettinessAK8:tau3"));
	//cout<<"set n subjettiness"<<endl;
        for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
            CAWDaughterPt     . push_back(ijet->daughter(ui)->pt());
            CAWDaughterEta    . push_back(ijet->daughter(ui)->eta());
            CAWDaughterPhi    . push_back(ijet->daughter(ui)->phi());
            CAWDaughterEnergy . push_back(ijet->daughter(ui)->energy());
            
            CAWDaughterMotherIndex . push_back(index);
        }
    }
    
    //Four std::vector
    SetValue("CAWJetPt"     , CAWJetPt);
    SetValue("CAWJetEta"    , CAWJetEta);
    SetValue("CAWJetPhi"    , CAWJetPhi);
    SetValue("CAWJetEnergy" , CAWJetEnergy);
    
    SetValue("CAWJetCSV"    , CAWJetCSV);
    //   SetValue("CAWJetRCN"    , CAWJetRCN);
    
    //Identity
    SetValue("CAWJetIndex"      , CAWJetIndex);
    SetValue("CAWJetnDaughters" , CAWJetnDaughters);
    
    //Mass
    SetValue("CAWJetTrimmedMass"     , CAWJetTrimmedMass);
    SetValue("CAWJetPrunedMass"     , CAWJetPrunedMass);
    SetValue("CAWJetFilteredMass"     , CAWJetFilteredMass);
    
    //Daughter four std::vector and index
    SetValue("CAWDaughterPt"     , CAWDaughterPt);
    SetValue("CAWDaughterEta"    , CAWDaughterEta);
    SetValue("CAWDaughterPhi"    , CAWDaughterPhi);
    SetValue("CAWDaughterEnergy" , CAWDaughterEnergy);
    
    SetValue("CAWDaughterMotherIndex" , CAWDaughterMotherIndex);
    
    //Get all CA8 jets (not just for W and Top)
    edm::InputTag CA8JetColl = edm::InputTag("slimmedJetsAK8");
    edm::Handle<std::vector<pat::Jet> > CA8Jets;
    event.getByLabel(CA8JetColl, CA8Jets);

    
    //Four std::vector
    std::vector <double> CA8JetPt;
    std::vector <double> CA8JetEta;
    std::vector <double> CA8JetPhi;
    std::vector <double> CA8JetEnergy;
    
    std::vector <double> CA8JetCSV;
    //   std::vector <double> CA8JetRCN;
    
    for (std::vector<pat::Jet>::const_iterator ijet = CA8Jets->begin(); ijet != CA8Jets->end(); ijet++){
        
        //Four std::vector
        CA8JetPt     . push_back(ijet->pt());
        CA8JetEta    . push_back(ijet->eta());
        CA8JetPhi    . push_back(ijet->phi());
        CA8JetEnergy . push_back(ijet->energy());
        
        CA8JetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
        //     CA8JetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));
    }
    
    //Four std::vector
    SetValue("CA8JetPt"     , CA8JetPt);
    SetValue("CA8JetEta"    , CA8JetEta);
    SetValue("CA8JetPhi"    , CA8JetPhi);
    SetValue("CA8JetEnergy" , CA8JetEnergy);
    
    SetValue("CA8JetCSV"    , CA8JetCSV);
    //   SetValue("CA8JetRCN"    , CA8JetRCN);
    
    //Get AK4 Jets
    //Four std::vector
    std::vector <double> AK4JetPt;
    std::vector <double> AK4JetEta;
    std::vector <double> AK4JetPhi;
    std::vector <double> AK4JetEnergy;
    
    std::vector <int>    AK4JetTBag;
    std::vector <double> AK4JetRCN;
    
    for (std::vector<edm::Ptr<pat::Jet> >::const_iterator ijet = vSelJets.begin();
         ijet != vSelJets.end(); ijet++){
        
        //Four std::vector
        TLorentzVector lv = selector->correctJet(**ijet, event);
        
        AK4JetPt     . push_back(lv.Pt());
        AK4JetEta    . push_back(lv.Eta());
        AK4JetPhi    . push_back(lv.Phi());
        AK4JetEnergy . push_back(lv.Energy());
        
        AK4JetTBag   . push_back(selector->isJetTagged(**ijet, event));
        AK4JetRCN    . push_back(((*ijet)->chargedEmEnergy()+(*ijet)->chargedHadronEnergy()) / ((*ijet)->neutralEmEnergy()+(*ijet)->neutralHadronEnergy()));
    }
    
    //Four std::vector
    SetValue("AK4JetPt"     , AK4JetPt);
    SetValue("AK4JetEta"    , AK4JetEta);
    SetValue("AK4JetPhi"    , AK4JetPhi);
    SetValue("AK4JetEnergy" , AK4JetEnergy);
    
    SetValue("AK4JetTBag"   , AK4JetTBag);
    SetValue("AK4JetRCN"    , AK4JetRCN);

    //Get cleaned AK4 Jets
    //Four std::vector
    std::vector <double> cleanedAK4JetPt;
    std::vector <double> cleanedAK4JetEta;
    std::vector <double> cleanedAK4JetPhi;
    std::vector <double> cleanedAK4JetEnergy;
    
    std::vector <int>    cleanedAK4JetTBag;
    std::vector <double> cleanedAK4JetRCN;
    
    for (std::vector<pat::Jet>::const_iterator ijet = vSelCleanedJets.begin();
         ijet != vSelCleanedJets.end(); ijet++){
      //no need to correct so just push back quantities from jet directly
        
      cleanedAK4JetPt     . push_back((*ijet).pt());
      cleanedAK4JetEta    . push_back((*ijet).eta());
      cleanedAK4JetPhi    . push_back((*ijet).phi());
      cleanedAK4JetEnergy . push_back((*ijet).energy());
      
      cleanedAK4JetTBag   . push_back(selector->isJetTagged(*ijet, event));
      cleanedAK4JetRCN    . push_back(((*ijet).chargedEmEnergy()+(*ijet).chargedHadronEnergy()) / ((*ijet).neutralEmEnergy()+(*ijet).neutralHadronEnergy()));
    }
    
    //Four std::vector
    SetValue("cleanedAK4JetPt"     , cleanedAK4JetPt);
    SetValue("cleanedAK4JetEta"    , cleanedAK4JetEta);
    SetValue("cleanedAK4JetPhi"    , cleanedAK4JetPhi);
    SetValue("cleanedAK4JetEnergy" , cleanedAK4JetEnergy);
    
    SetValue("cleanedAK4JetTBag"   , cleanedAK4JetTBag);
    SetValue("cleanedAK4JetRCN"    , cleanedAK4JetRCN);

    
    // MET
    double _met = -9999.0;
    double _met_phi = -9999.0;
    // Corrected MET
    double _corr_met = -9999.0;
    double _corr_met_phi = -9999.0;
    
    if(pMet.isNonnull() && pMet.isAvailable()) {
        _met = pMet->p4().pt();
        _met_phi = pMet->p4().phi();
        
        TLorentzVector corrMET = selector->correctMet(*pMet, event, vSelCleanedJets);
        if(corrMET.Pt()>0) {
            _corr_met = corrMET.Pt();
            _corr_met_phi = corrMET.Phi();
        }
        
    }
    SetValue("met", _met);
    SetValue("met_phi", _met_phi);
    SetValue("corr_met", _corr_met);
    SetValue("corr_met_phi", _corr_met_phi);
    
    
    //
    //_____ Gen Info ______________________________
    //
    
    //Four std::vector
    std::vector <double> genPt;
    std::vector <double> genEta;
    std::vector <double> genPhi;
    std::vector <double> genEnergy;
    
    //Identity
    std::vector <int> genID;
    std::vector <int> genIndex;
    std::vector <int> genStatus;
    std::vector <int> genMotherID;
    std::vector <int> genMotherIndex;
    std::vector<bool> genIsFromTau;
    std::vector<bool> genIsPrompt;
    std::vector<bool> genPMotherHasC;
    std::vector<bool> genPMotherHasB;
    std::vector<int> genPMother;
    //event weights
    std::vector<double> evtWeightsMC;
    float MCWeight=1;
    std::vector<double> LHEweights;
    std::vector<int> LHEweightids;

    if (isMc){
      //load info for event weight
      edm::Handle<GenEventInfoProduct> genEvtInfo;
      edm::InputTag gen_it("generator");
      event.getByLabel(gen_it, genEvtInfo );

      std::vector<double> evtWeights = genEvtInfo->weights();
      double theWeight = genEvtInfo->weight();
      
      evtWeightsMC=evtWeights;
      MCWeight = theWeight;
   

      //weights for mc uncertainties
      edm::Handle<LHEEventProduct> EvtHandle;
      edm::InputTag theSrc("externalLHEProducer");
      if(event.getByLabel(theSrc,EvtHandle)){
	  
	// Storing LHE weights https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
	// for MC@NLO renormalization and factorization scale. 
	// ID numbers 1001 - 1009. (muR,muF) = 
	// 0 = 1001: (1,1)    3 = 1004: (2,1)    6 = 1007: (0.5,1)  
	// 1 = 1002: (1,2)    4 = 1005: (2,2)   7 = 1008: (0.5,2)  
	// 2 = 1003: (1,0.5)  5 = 1006: (2,0.5) 8 = 1009: (0.5,0.5)
	// for PDF variations: ID numbers > 2000

	std::string weightidstr;
	int weightid;
	if(EvtHandle->weights().size() > 0){  
	  for(unsigned int i = 0; i < EvtHandle->weights().size(); i++){
	    weightidstr = EvtHandle->weights()[i].id;
	    weightid = std::stoi(weightidstr);
	    LHEweights.push_back(EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP());
	    LHEweightids.push_back(weightid);
	  }
	}
      }

      //load genparticles collection
      edm::Handle<reco::GenParticleCollection> genParticles;
      event.getByLabel(genParticles_it, genParticles);
        
        for(size_t i = 0; i < genParticles->size(); i++){
            const reco::GenParticle & p = (*genParticles).at(i);

	    bool keep=false;

	    //attempt using mc gen particle flags
	    if(p.status()==1) keep = true; //all stable out going particles
	    else if( p.isPromptDecayed() ) keep = true; //keep prompt particles which have undergone decay (e.g. B-mesons)
	    else if( p.statusFlags().isPrompt()) keep = true; // keep all prompt particles in case I want to do quark matching
	    else if( p.isDirectPromptTauDecayProductFinalState()) keep = true; //log save leptons from tau decays since they are 'prompt' for us
	    

	    /* OLD METHOD USING STATUSES IS FAULTY
            //Find particles from hard scattering/stable
            if (p.status() >=20 && p.status()<=29) keep=true;
	    else if(p.status()==1 && ( fabs(p.pdgId())>10 && fabs(p.pdgId())<19) ){
	      //since we know it's a status 1 lepton, only keep if it's mother is a W and either status 22 or 52
	      if(p.mother()){
		reco::Candidate* mother = (reco::Candidate*) p.mother();
		if( fabs(mother->pdgId())==24 && ( (mother->status()>=20 && mother->status()<=29) || (mother->status()>=50 && mother->status()<=59)) ) keep=true;
	      }
	    }*/



            if(!keep) continue;    

	    bool promptMotherHasB = false;
	    bool promptMotherHasC = false;

	    reco::GenParticle* mother = 0;
	    if(p.status()==1){
	      mother = (reco::GenParticle*) p.mother();
	      while(mother){
		if(mother->isPromptDecayed()) break;
		else{ mother = (reco::GenParticle*) mother->mother();}
	      }
	    }

	    if(mother){
	      HepPDT::ParticleID heppdtid(mother->pdgId());
	      promptMotherHasB = heppdtid.hasBottom();
	      promptMotherHasC = heppdtid.hasCharm();
	    }

	    genPMotherHasB.push_back(promptMotherHasB);
	    genPMotherHasC.push_back(promptMotherHasC);
	    if(mother) genPMother.push_back(mother->pdgId());
	    else       genPMother.push_back(-999);
	    //Four std::vector
	    genPt     . push_back(p.pt());
	    genEta    . push_back(p.eta());
	    genPhi    . push_back(p.phi());
	    genEnergy . push_back(p.energy());
            
	    //Identity
	    genID            . push_back(p.pdgId());
	    genIndex         . push_back((int) i);
	    genStatus        . push_back(p.status());
	    genIsPrompt      . push_back(p.statusFlags().isPrompt());
	    genIsFromTau     . push_back(p.isDirectPromptTauDecayProductFinalState()); //save whether or not a particle came from a prompt tau
	    //mother info - stores mother pdgID if mother exists or -999 if mother doesn't exist. This way the entries will always be aligned with main gen particle entries
	    if(p.mother()) genMotherID.push_back((p.mother())->pdgId());
	    else genMotherID.push_back(-999);
	    //genMotherIndex   . push_back(mInd);
	
	}//End loop over gen particles
    }  //End MC-only if
    
    // Four std::vector
    SetValue("genPt"    , genPt);
    SetValue("genEta"   , genEta);
    SetValue("genPhi"   , genPhi);
    SetValue("genEnergy", genEnergy);
    
    // Identity
    SetValue("genID"         , genID);
    SetValue("genIndex"      , genIndex);
    SetValue("genStatus"     , genStatus);
    SetValue("genMotherID"   , genMotherID);
    SetValue("genMotherIndex", genMotherIndex);
    SetValue("genIsPrompt"   , genIsPrompt);
    SetValue("genIsFromPromptTau", genIsFromTau);
    SetValue("evtWeightsMC", evtWeightsMC);
    SetValue("MCWeight",MCWeight);
    SetValue("LHEWeights",LHEweights);
    SetValue("LHEWeightIDs",LHEweightids);
    SetValue("genPMotherHasC",genPMotherHasC);
    SetValue("genPMotherHasB",genPMotherHasB);
    SetValue("genPMother",genPMother);

    return 0;
}

int DileptonCalc::findMatch(const reco::GenParticleCollection & genParticles, int idToMatch, double eta, double phi)
{
    float dRtmp = 1000;
    float closestDR = 10000.;
    int closestGenPart = -1;
    
    for(size_t j = 0; j < genParticles.size(); ++ j) {
        const reco::GenParticle & p = (genParticles).at(j);
        dRtmp = mdeltaR( eta, phi, p.eta(), p.phi());
        if ( dRtmp < closestDR && abs(p.pdgId()) == idToMatch){// && dRtmp < 0.3) {
            closestDR = dRtmp;
            closestGenPart = j;
        }//end of requirement for matching
    }//end of gen particle loop 
    return closestGenPart;
}


double DileptonCalc::mdeltaR(double eta1, double phi1, double eta2, double phi2) {
    return std::sqrt(deltaR2 (eta1, phi1, eta2, phi2));
}

void DileptonCalc::fillMotherInfo(const reco::Candidate *mother, int i, std::vector <int> & momid, std::vector <int> & momstatus, std::vector<double> & mompt, std::vector<double> & mometa, std::vector<double> & momphi, std::vector<double> & momenergy)
{
    if(mother) {
        momid.push_back(mother->pdgId());
        momstatus.push_back(mother->status());
        mompt.push_back(mother->pt());
        mometa.push_back(mother->eta());
        momphi.push_back(mother->phi());
        momenergy.push_back(mother->energy());
        if(i<10)fillMotherInfo(mother->mother(), i+1, momid, momstatus, mompt, mometa, momphi, momenergy);
    }
    
    
} 
