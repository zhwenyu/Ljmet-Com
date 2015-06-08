/*
 Calculator for substructure variables
 
 Author: Joshua Swanson, 2014. Updated Julie Hogan, 2015. 
 */

#include <iostream>
#include <limits>   // std::numeric_limits
#include <vector>
#include <string>
#include "TLorentzVector.h"

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"
#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h"

class LjmetFactory;

class JetSubCalc : public BaseCalc {
    
public:
    JetSubCalc();
    virtual ~JetSubCalc();
    virtual int BeginJob();
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob();
    
private:
  edm::InputTag slimmedJetColl_it;
  edm::InputTag slimmedJetsAK8Coll_it;
  edm::InputTag selectedPatJetsAK8Coll_it;
  edm::InputTag selectedPatSubJetsAK8Coll_it;
  edm::InputTag selectedPatPackedSubJetsAK8Coll_it;
  edm::InputTag selectedPatJetsCA15Coll_it;
  edm::InputTag httTagInfo_it;
  std::string bDiscriminant;
  std::string tagInfo;
};

static int reg = LjmetFactory::GetInstance()->Register(new JetSubCalc(), "JetSubCalc");

JetSubCalc::JetSubCalc()
{
}

JetSubCalc::~JetSubCalc()
{
}

int JetSubCalc::BeginJob()
{
    if (mPset.exists("slimmedJetColl")) slimmedJetColl_it = mPset.getParameter<edm::InputTag>("slimmedJetColl");
    else slimmedJetColl_it = edm::InputTag("slimmedJets");
    
    if (mPset.exists("slimmedJetsAK8Coll")) slimmedJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("slimmedJetsAK8Coll");
    else slimmedJetsAK8Coll_it = edm::InputTag("slimmedJetsAK8");
    
    if (mPset.exists("bDiscriminant")) bDiscriminant = mPset.getParameter<std::string>("bDiscriminant");
    else bDiscriminant = "combinedInclusiveSecondaryVertexV2BJetTags";
    
    std::cout << " JetSubCalc Bdisc = " << bDiscriminant << std::endl;

    if (mPset.exists("tagInfo")) tagInfo = mPset.getParameter<std::string>("tagInfo");
    else tagInfo = "caTop";

    // Comment out these inclusions for files without HTT, SoftDrop, Subjets, etc.
    if (mPset.exists("httTagInfo")) httTagInfo_it = mPset.getParameter<edm::InputTag>("httTagInfo");
    else httTagInfo_it = edm::InputTag("HTT");

    if (mPset.exists("selectedJetsAK8Coll")) selectedPatJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("selectedJetsAK8Coll");
    else selectedPatJetsAK8Coll_it = edm::InputTag("selectedPatJetsAK8PFCHS");

    if (mPset.exists("selectedJetsCA15Coll")) selectedPatJetsCA15Coll_it = mPset.getParameter<edm::InputTag>("selectedJetsCA15Coll");
    else selectedPatJetsCA15Coll_it = edm::InputTag("selectedPatJetsCA15PFCHS");

    if (mPset.exists("selectedSubJetsAK8Coll")) selectedPatSubJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("selectedSubJetsAK8Coll");
    else selectedPatSubJetsAK8Coll_it = edm::InputTag("selectedPatJetsAK8PFCHSSoftDropSubjets");

    if (mPset.exists("selectedPackedSubJetsAK8Coll")) selectedPatPackedSubJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("selectedPackedSubJetsAK8Coll");
    else selectedPatPackedSubJetsAK8Coll_it = edm::InputTag("selectedPatJetsAK8PFCHSSoftDropPacked");
    
    return 0;
}

int JetSubCalc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector)
{
    // ----- Get AK4 jet objects from the selector -----
    // This is updated -- original version used all AK4 jets without selection 
    std::vector<edm::Ptr<pat::Jet> >            const & theJets = selector->GetSelectedJets();
    std::vector<std::pair<TLorentzVector,bool>> const & theCorrBtagJets = selector->GetCorrJetsWithBTags();

    // Old version
    //    edm::Handle<std::vector<pat::Jet> > theJets;
    //    event.getByLabel(slimmedJetColl_it, theJets);

    float subjetCSV;
    int CSVL, CSVM, CSVT;
    double theJetHT;
    
    // Available variables
    std::vector<double> theJetPt;
    std::vector<double> theJetEta;
    std::vector<double> theJetPhi;
    std::vector<double> theJetEnergy;
    std::vector<double> theJetCSV;
    
    // Additional variables related to the associated secondary vertex if there is one
    // Mass of the vertex
    std::vector<double> theJetVtxMass;
    // Number of tracks
    std::vector<double> theJetVtxNtracks;
    // Decay length value and significance
    std::vector<double> theJetVtx3DVal;
    std::vector<double> theJetVtx3DSig;
    
    // Discriminator for the MVA PileUp id.
    // NOTE: Training used is for ak5PFJetsCHS in CMSSW 5.3.X and Run 1 pileup
    std::vector<double> theJetPileupJetId;
    
    //Identity
    std::vector<int> theJetIndex;
    std::vector<int> theJetnDaughters;
    
    //Daughter four vector and index
    std::vector<double> theJetDaughterPt;
    std::vector<double> theJetDaughterEta;
    std::vector<double> theJetDaughterPhi;
    std::vector<double> theJetDaughterEnergy;
    
    std::vector<int> theJetDaughterMotherIndex;
    
    std::vector<int> theJetCSVLSubJets;
    std::vector<int> theJetCSVMSubJets;
    std::vector<int> theJetCSVTSubJets;

    std::vector<int> theJetFlav;
    std::vector<int> theJetBTag;
    
    double theVtxMass, theVtxNtracks, theVtx3DVal, theVtx3DSig, thePileupJetId;

    for (unsigned int ii = 0; ii < theCorrBtagJets.size(); ii++){
      int index = ii;
      
      theVtxNtracks  = -std::numeric_limits<double>::max();
      theVtxNtracks  = (double)theJets[ii]->userFloat("vtxNtracks");
				
      if (theVtxNtracks > 0) {
	theVtxMass     = -std::numeric_limits<double>::max();
	theVtx3DVal    = -std::numeric_limits<double>::max();
	theVtx3DSig    = -std::numeric_limits<double>::max();
	
	theVtxMass     = (double)theJets[ii]->userFloat("vtxMass");
	theVtx3DVal    = (double)theJets[ii]->userFloat("vtx3DVal");
	theVtx3DSig    = (double)theJets[ii]->userFloat("vtx3DSig");
        
	theJetVtxNtracks.push_back(theVtxNtracks);
	theJetVtxMass.push_back(theVtxMass);
	theJetVtx3DVal.push_back(theVtx3DVal);
	theJetVtx3DSig.push_back(theVtx3DSig);
      }
      
      thePileupJetId = -std::numeric_limits<double>::max();
      thePileupJetId = (double)theJets[ii]->userFloat("pileupJetId:fullDiscriminant");
      theJetPileupJetId.push_back(thePileupJetId);

      //Four vector
      TLorentzVector lv = theCorrBtagJets[ii].first;

      theJetPt     . push_back(lv.Pt());
      theJetEta    . push_back(lv.Eta());
      theJetPhi    . push_back(lv.Phi());
      theJetEnergy . push_back(lv.Energy());
      
      theJetCSV.push_back(theJets[ii]->bDiscriminator( bDiscriminant ));
      theJetBTag.push_back(theCorrBtagJets[ii].second);
      theJetFlav.push_back(abs(theJets[ii]->partonFlavour()));
      
      theJetIndex.push_back(index);
      theJetnDaughters.push_back((int)theJets[ii]->numberOfDaughters());
 
      //HT
      theJetHT += lv.Pt(); 
    }

    /* // Old version
    for (std::vector<edm::Ptr<pat::Jet> >::const_iterator ijet = theJets.begin(); ijet != theJets.end(); ijet++) {
      int index = (int)((*ijet)-theJets.begin());
      
      theVtxNtracks  = -std::numeric_limits<double>::max();
      theVtxNtracks  = (double)(*ijet)->userFloat("vtxNtracks");
      
      if (theVtxNtracks > 0) {
	theVtxMass     = -std::numeric_limits<double>::max();
	theVtx3DVal    = -std::numeric_limits<double>::max();
	theVtx3DSig    = -std::numeric_limits<double>::max();
	
	theVtxMass     = (double)(*ijet)->userFloat("vtxMass");
	theVtx3DVal    = (double)(*ijet)->userFloat("vtx3DVal");
	theVtx3DSig    = (double)(*ijet)->userFloat("vtx3DSig");
        
	theJetVtxNtracks.push_back(theVtxNtracks);
	theJetVtxMass.push_back(theVtxMass);
	theJetVtx3DVal.push_back(theVtx3DVal);
	theJetVtx3DSig.push_back(theVtx3DSig);
      }
      
      thePileupJetId = -std::numeric_limits<double>::max();
      thePileupJetId = (double)(*ijet)->userFloat("pileupJetId:fullDiscriminant");
      theJetPileupJetId.push_back(thePileupJetId);
      
      theJetPt    .push_back((*ijet)->pt());
      theJetEta   .push_back((*ijet)->eta());
      theJetPhi   .push_back((*ijet)->phi());
      theJetEnergy.push_back((*ijet)->energy());
      theJetCSV.push_back((*ijet)->bDiscriminator( bDiscriminant ));
      
      theJetIndex.push_back(index);
      theJetnDaughters.push_back((int)(*ijet)->numberOfDaughters());
      
      CSVL = 0;
      CSVM = 0;
      CSVT = 0;
      subjetCSV = -std::numeric_limits<float>::max();
      
      for (size_t ui = 0; ui < (*ijet)->numberOfDaughters(); ui++) {
	pat::PackedCandidate const * theDaughter = dynamic_cast<pat::PackedCandidate const *>((*ijet)->daughter(ui));
        
	theJetDaughterPt    .push_back(theDaughter->pt());
	theJetDaughterEta   .push_back(theDaughter->eta());
	theJetDaughterPhi   .push_back(theDaughter->phi());
	theJetDaughterEnergy.push_back(theDaughter->energy());
        
	theJetDaughterMotherIndex.push_back(index);
        
        
	//            pat::Jet const * subjet = dynamic_cast<pat::Jet const *>((*ijet)->daughter(ui));
	//            subjetCSV = subjet->bDiscriminator(bDiscriminant);
	//            if (subjetCSV > 0.244 && subjet->pt() > 20.) {
	//                CSVL++;
	//            }
	//            if (subjetCSV > 0.679 && subjet->pt() > 20.) {
	//                CSVM++;
	//            }
	//            if (subjetCSV > 0.898 && subjet->pt() > 20.) {
	//                CSVT++;
	//            }
      }
      theJetCSVLSubJets.push_back(CSVL);
      theJetCSVMSubJets.push_back(CSVM);
      theJetCSVTSubJets.push_back(CSVT);
    }
    */
    SetValue("theJetPt",     theJetPt);
    SetValue("theJetEta",    theJetEta);
    SetValue("theJetPhi",    theJetPhi);
    SetValue("theJetEnergy", theJetEnergy);
    SetValue("theJetCSV",    theJetCSV);
    SetValue("theJetBTag",   theJetBTag);
    SetValue("theJetFlav",   theJetFlav);
    
    SetValue("theJetVtxMass",     theJetVtxMass);
    SetValue("theJetVtxNtracks",  theJetVtxNtracks);
    SetValue("theJetVtx3DVal",    theJetVtx3DVal);
    SetValue("theJetVtx3DSig",    theJetVtx3DSig);
    SetValue("theJetPileupJetId", theJetPileupJetId);
    
    //    SetValue("theJetIndex",      theJetIndex);
    SetValue("theJetnDaughters", theJetnDaughters);
    /*    
    SetValue("theJetDaughterPt",     theJetDaughterPt);
    SetValue("theJetDaughterEta",    theJetDaughterEta);
    SetValue("theJetDaughterPhi",    theJetDaughterPhi);
    SetValue("theJetDaughterEnergy", theJetDaughterEnergy);
    */
    //    SetValue("theJetDaughterMotherIndex", theJetDaughterMotherIndex);
    
    //    SetValue("theJetCSVLSubJets", theJetCSVLSubJets);
    //    SetValue("theJetCSVMSubJets", theJetCSVMSubJets);
    //    SetValue("theJetCSVTSubJets", theJetCSVTSubJets);
    
    // Load in AK8 jets (no selection performed on these)
    edm::Handle<std::vector<pat::Jet> > theAK8Jets;
    event.getByLabel(slimmedJetsAK8Coll_it, theAK8Jets);

    // Uncomment these for standard miniAOD
    edm::Handle<std::vector<pat::Jet> > theAK8JetsSD;
    event.getByLabel(selectedPatJetsAK8Coll_it, theAK8JetsSD);

    edm::Handle<std::vector<pat::Jet> > theAK8JetsSDsubjets;
    event.getByLabel(selectedPatSubJetsAK8Coll_it, theAK8JetsSDsubjets);

    edm::Handle<std::vector<pat::Jet> > theAK8JetsSDPackedsubjets;
    event.getByLabel(selectedPatPackedSubJetsAK8Coll_it, theAK8JetsSDPackedsubjets);
    
    // Four vector
    std::vector<double> theJetAK8Pt;
    std::vector<double> theJetAK8Eta;
    std::vector<double> theJetAK8Phi;
    std::vector<double> theJetAK8Energy;
    std::vector<double> theJetAK8CSV;
    
    // Pruned, trimmed and filtered masses available
    std::vector<double> theJetAK8PrunedMass;
    std::vector<double> theJetAK8TrimmedMass;
    std::vector<double> theJetAK8FilteredMass;
    std::vector<double> theJetAK8SoftDropMass;
    
    // n-subjettiness variables tau1, tau2, and tau3 available
    std::vector<double> theJetAK8NjettinessTau1;
    std::vector<double> theJetAK8NjettinessTau2;
    std::vector<double> theJetAK8NjettinessTau3;
    
    std::vector<double> theJetAK8caTopTopMass;
    std::vector<double> theJetAK8caTopMinMass;
    std::vector<int> theJetAK8caTopnSubJets;

    std::vector<double> theJetAK8Mass;
    std::vector<int>    theJetAK8Index;
    std::vector<int>    theJetAK8nDaughters;
    
    // Daughter four vector and index
    std::vector<double> theJetAK8DaughterPt;
    std::vector<double> theJetAK8DaughterEta;
    std::vector<double> theJetAK8DaughterPhi;
    std::vector<double> theJetAK8DaughterEnergy;
    
    std::vector<int> theJetAK8DaughterMotherIndex;
    
    std::vector<int> theJetAK8CSVLSubJets;
    std::vector<int> theJetAK8CSVMSubJets;
    std::vector<int> theJetAK8CSVTSubJets;

    std::vector<int> theJetAK8SoftDropNsubjets;
    std::vector<double> theJetAK8SoftDropSubjet1CSV;
    std::vector<double> theJetAK8SoftDropSubjet2CSV;
    
    double topMass, minMass, theHTTTopMass, theHTTfrec;
    int nSubJets, numbersubjets;
    double thePrunedMass, theTrimmedMass, theFilteredMass, theSoftDropMass;
    double theNjettinessTau1, theNjettinessTau2, theNjettinessTau3;
    double sub1csv, sub2csv;

    // Check sizes -- comment out for standard miniAOD
    if(theAK8Jets->size() != theAK8JetsSD->size()){
      std::cout << "Jet collection size mismatch: slimmed = " << theAK8Jets->size() << ", softdrop = " << theAK8JetsSD->size() << std::endl;
    }

    for (std::vector<pat::Jet>::const_iterator ijet = theAK8Jets->begin(); ijet != theAK8Jets->end(); ijet++) {
      int index = (int)(ijet-theAK8Jets->begin());
      
      thePrunedMass   = -std::numeric_limits<double>::max();
      theTrimmedMass  = -std::numeric_limits<double>::max();
      theFilteredMass = -std::numeric_limits<double>::max();
      
      thePrunedMass   = (double)ijet->userFloat("ak8PFJetsCHSPrunedLinks");
      theTrimmedMass  = (double)ijet->userFloat("ak8PFJetsCHSTrimmedLinks");
      theFilteredMass = (double)ijet->userFloat("ak8PFJetsCHSFilteredLinks");
      
      theNjettinessTau1 = -std::numeric_limits<double>::max();
      theNjettinessTau2 = -std::numeric_limits<double>::max();
      theNjettinessTau3 = -std::numeric_limits<double>::max();
      theNjettinessTau1 = (double)ijet->userFloat("NjettinessAK8:tau1");
      theNjettinessTau2 = (double)ijet->userFloat("NjettinessAK8:tau2");
      theNjettinessTau3 = (double)ijet->userFloat("NjettinessAK8:tau3");
      
      theJetAK8Pt    .push_back(ijet->pt());
      theJetAK8Eta   .push_back(ijet->eta());
      theJetAK8Phi   .push_back(ijet->phi());
      theJetAK8Energy.push_back(ijet->energy());
      theJetAK8CSV.push_back(ijet->bDiscriminator( bDiscriminant ));
      
      theJetAK8PrunedMass  .push_back(thePrunedMass);
      theJetAK8TrimmedMass .push_back(theTrimmedMass);
      theJetAK8FilteredMass.push_back(theFilteredMass);
      
      theJetAK8NjettinessTau1.push_back(theNjettinessTau1);
      theJetAK8NjettinessTau2.push_back(theNjettinessTau2);
      theJetAK8NjettinessTau3.push_back(theNjettinessTau3);
      
      theJetAK8Mass.push_back(ijet->mass());
      theJetAK8nDaughters.push_back((int)ijet->numberOfDaughters());
      
      theJetAK8Index.push_back(index);
      
      reco::CATopJetTagInfo const * jetInfo = dynamic_cast<reco::CATopJetTagInfo const *>( ijet->tagInfo( tagInfo ));
      
      if ( jetInfo != 0 ) {
	topMass   = -std::numeric_limits<double>::max();
	minMass   = -std::numeric_limits<double>::max();
	nSubJets  = std::numeric_limits<int>::min();
        
	topMass = jetInfo->properties().topMass;
	minMass = jetInfo->properties().minMass;
	nSubJets = jetInfo->properties().nSubJets;
        
	theJetAK8caTopTopMass.push_back(topMass);
	theJetAK8caTopMinMass.push_back(minMass);
	theJetAK8caTopnSubJets.push_back(nSubJets);
      }
      
      subjetCSV = -std::numeric_limits<float>::max();
      CSVL = 0;
      CSVM = 0;
      CSVT = 0;
      
      for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++) {
	pat::PackedCandidate const * theDaughter = dynamic_cast<pat::PackedCandidate const *>(ijet->daughter(ui));
	theJetAK8DaughterPt    .push_back(theDaughter->pt());
	theJetAK8DaughterEta   .push_back(theDaughter->eta());
	theJetAK8DaughterPhi   .push_back(theDaughter->phi());
	theJetAK8DaughterEnergy.push_back(theDaughter->energy());
        
	theJetAK8DaughterMotherIndex.push_back(index);
	/*            
		      subjetCSV = subjet->bDiscriminator(bDiscriminant);
		      if (subjetCSV > 0.244 && ijet->daughter(ui)->pt() > 20.) {
		      CSVL++;
		      }
		      if (subjetCSV > 0.679 && ijet->daughter(ui)->pt() > 20.){
		      CSVM++;
		      }
		      if (subjetCSV > 0.898 && ijet->daughter(ui)->pt() > 20.){
		      CSVT++;
		      }
	*/
      }
      theJetAK8CSVLSubJets.push_back(CSVL);
      theJetAK8CSVMSubJets.push_back(CSVM);
      theJetAK8CSVTSubJets.push_back(CSVT);

      // Match the slimmedAK8 jet to the selectedPatAK8 jet (only check DR if collection sizes don't match)
      // Comment out for standard miniAOD
      std::vector<pat::Jet>::const_iterator ijetSD = theAK8JetsSD->begin()+index;
      if(theAK8Jets->size() != theAK8JetsSD->size()){
	double minDR = 100;
	for (std::vector<pat::Jet>::const_iterator ijet2 = theAK8JetsSD->begin(); ijet2 != theAK8JetsSD->end(); ijet2++) {
	  TLorentzVector jet1; TLorentzVector jet2;
	  jet1.SetPtEtaPhiE(ijet->pt(),ijet->eta(),ijet->phi(),ijet->energy());
	  jet2.SetPtEtaPhiE(ijet2->pt(),ijet2->eta(),ijet2->phi(),ijet2->energy());
	  double DR = jet1.DeltaR(jet2);
	  if(DR < minDR){
	    minDR = DR;
	    ijetSD = ijet2;
	  }
	}
	if(minDR > 0.8) std::cout << "Didn't find match within cone size 0.8: minDR = " << minDR << std::endl;
      }

      theSoftDropMass = -std::numeric_limits<double>::max();
      if(ijetSD != theAK8JetsSD->end()){
	theSoftDropMass = (double)ijetSD->userFloat("ak8PFJetsCHSSoftDropLinks");
	theJetAK8SoftDropMass.push_back(theSoftDropMass);
      }
      // Match the selectedPatAK8 jet to a collection of packed subjets using DR
      // comment out for standard miniAOD
      double minDRpack = 1000;
      std::vector<pat::Jet>::const_iterator ijetPacked = theAK8JetsSDPackedsubjets->begin();
      for (std::vector<pat::Jet>::const_iterator pack = theAK8JetsSDPackedsubjets->begin(); pack != theAK8JetsSDPackedsubjets->end(); pack++) {

	TLorentzVector packedsjs; TLorentzVector sdjet;
	sdjet.SetPtEtaPhiE(ijetSD->pt(),ijetSD->eta(),ijetSD->phi(),ijetSD->energy());
	packedsjs.SetPtEtaPhiE(pack->pt(),pack->eta(),pack->phi(),pack->energy());

	double DRpack = sdjet.DeltaR(packedsjs);
	if(DRpack < minDRpack){
	  minDRpack = DRpack;
	  ijetPacked = pack;
	}
      }

      numbersubjets = -std::numeric_limits<int>::max();
      if(minDRpack > 0.8){
	//	std::cout << "Didn't find packed subjets within cone size 0.8" << std::endl;
	numbersubjets = 0;
      }else{ numbersubjets = ijetPacked->numberOfDaughters(); }
      theJetAK8SoftDropNsubjets.push_back(numbersubjets);
      
      // If there are at least 2 subjets, get first 2 vectors and match to the subjet collection for B-discriminant
      std::vector<pat::Jet>::const_iterator subjet1 = theAK8JetsSDsubjets->begin();
      std::vector<pat::Jet>::const_iterator subjet2 = theAK8JetsSDsubjets->begin();
      sub1csv = -std::numeric_limits<double>::max();
      sub2csv = -std::numeric_limits<double>::max();

      int gotsub[2] = {1,1};
      double minDRsubjet = 1000;
      if(numbersubjets > 1){
	for (int ui = 0; ui < 2; ui++) {
	  auto const * theDaughter = ijetPacked->daughter(ui);

	  TLorentzVector daughter;
	  daughter.SetPtEtaPhiE(theDaughter->pt(),theDaughter->eta(),theDaughter->phi(),theDaughter->energy());

	  minDRsubjet = 1000;
	  for (std::vector<pat::Jet>::const_iterator sub = theAK8JetsSDsubjets->begin(); sub != theAK8JetsSDsubjets->end(); sub++) {
	    TLorentzVector sj;
	    sj.SetPtEtaPhiE(sub->pt(),sub->eta(),sub->phi(),sub->energy());
	    double DR = daughter.DeltaR(sj);

	    if(DR < minDRsubjet){
	      minDRsubjet = DR;
	      if(ui == 0) subjet1 = sub;
	      if(ui == 1) subjet2 = sub;
	    }
	  }
	  if(minDRsubjet > 0.8){
	    std::cout << "Didn't find subjet match for Packed daughter " << ui << std::endl;
	    gotsub[ui] = 0;
	  }
	}

	if(gotsub[0] > 0) sub1csv = subjet1->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	else sub1csv = 0;

	if(gotsub[1] > 0) sub2csv = subjet2->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	else sub2csv = 0;
      }

      theJetAK8SoftDropSubjet1CSV.push_back(sub1csv);
      theJetAK8SoftDropSubjet2CSV.push_back(sub2csv);

      // Comment out all the way to here for standard miniAOD

    }

    SetValue("theJetAK8Pt",     theJetAK8Pt);
    SetValue("theJetAK8Eta",    theJetAK8Eta);
    SetValue("theJetAK8Phi",    theJetAK8Phi);
    SetValue("theJetAK8Energy", theJetAK8Energy);
    SetValue("theJetAK8CSV",    theJetAK8CSV);
    
    SetValue("theJetAK8PrunedMass",   theJetAK8PrunedMass);
    SetValue("theJetAK8TrimmedMass",  theJetAK8TrimmedMass);
    SetValue("theJetAK8FilteredMass", theJetAK8FilteredMass);
    SetValue("theJetAK8SoftDropMass", theJetAK8SoftDropMass);
    
    SetValue("theJetAK8NjettinessTau1", theJetAK8NjettinessTau1);
    SetValue("theJetAK8NjettinessTau2", theJetAK8NjettinessTau2);
    SetValue("theJetAK8NjettinessTau3", theJetAK8NjettinessTau3);
    
    SetValue("theJetAK8Mass",   theJetAK8Mass);
    
    //    SetValue("theJetAK8Index",      theJetAK8Index);
    SetValue("theJetAK8nDaughters", theJetAK8nDaughters);
    
    SetValue("theJetAK8caTopTopMass", theJetAK8caTopTopMass);
    SetValue("theJetAK8caTopMinMass", theJetAK8caTopMinMass);
    SetValue("theJetAK8caTopnSubJets", theJetAK8caTopnSubJets);
    /*    
    SetValue("theJetAK8DaughterPt",     theJetAK8DaughterPt);
    SetValue("theJetAK8DaughterEta",    theJetAK8DaughterEta);
    SetValue("theJetAK8DaughterPhi",    theJetAK8DaughterPhi);
    SetValue("theJetAK8DaughterEnergy", theJetAK8DaughterEnergy);
    */
    //    SetValue("theJetAK8DaughterMotherIndex", theJetAK8DaughterMotherIndex);
    
    //    SetValue("theJetAK8CSVLSubJets", theJetAK8CSVLSubJets);
    //    SetValue("theJetAK8CSVMSubJets", theJetAK8CSVMSubJets);
    //    SetValue("theJetAK8CSVTSubJets", theJetAK8CSVTSubJets);

    SetValue("theJetAK8SoftDropNsubjets", theJetAK8SoftDropNsubjets);
    SetValue("theJetAK8SoftDropSubjet1CSV", theJetAK8SoftDropSubjet1CSV);
    SetValue("theJetAK8SoftDropSubjet2CSV", theJetAK8SoftDropSubjet2CSV);

    // Comment out entire CA15 loop for standard miniAOD

    std::vector<double> theJetCA15Pt;
    std::vector<double> theJetCA15Eta;
    std::vector<double> theJetCA15Phi;
    std::vector<double> theJetCA15Energy;

    std::vector<double> theJetCA15HTTTopMass;
    std::vector<double> theJetCA15HTTfrec;
    std::vector<double> theJetCA15NjettinessTau1;
    std::vector<double> theJetCA15NjettinessTau2;
    std::vector<double> theJetCA15NjettinessTau3;

    edm::Handle<std::vector<pat::Jet> > theCA15Jets;
    event.getByLabel(selectedPatJetsCA15Coll_it, theCA15Jets);

    edm::Handle<std::vector<reco::HTTTopJetTagInfo> > httTagInfo;
    event.getByLabel(httTagInfo_it, httTagInfo);

    for (std::vector<pat::Jet>::const_iterator ijet = theCA15Jets->begin(); ijet != theCA15Jets->end(); ijet++) {

      theJetCA15Pt.push_back(ijet->pt());
      theJetCA15Eta.push_back(ijet->eta());
      theJetCA15Phi.push_back(ijet->phi());
      theJetCA15Energy.push_back(ijet->energy());

      TLorentzVector ca15;
      ca15.SetPtEtaPhiE(ijet->pt(),ijet->eta(),ijet->phi(),ijet->energy());

      double minDR = 1000;
      reco::HTTTopJetProperties CA15HTTProperties;
      for (std::vector<reco::HTTTopJetTagInfo>::const_iterator ijet = httTagInfo->begin(); ijet != httTagInfo->end(); ijet++) {
      
	reco::HTTTopJetTagInfo const jetInfo = (*ijet);
	reco::HTTTopJetProperties HTTProperties = jetInfo.properties();

	TLorentzVector fj;
	fj.SetPtEtaPhiM(HTTProperties.fjPt, HTTProperties.fjEta, HTTProperties.fjPhi, HTTProperties.fjMass);

	double DR = ca15.DeltaR(fj);
	if(DR < minDR){
	  minDR = DR;
	  CA15HTTProperties = HTTProperties;
	}
      }

      theHTTTopMass = -std::numeric_limits<double>::max();
      theHTTfrec = -std::numeric_limits<double>::max();
      if(minDR < 1.5){
	theHTTTopMass = CA15HTTProperties.topMass;
	theHTTfrec = CA15HTTProperties.fRec;
	//	std::cout << "Matched the CA15 jet to an HTT input fatjet with cone size 1.5" << std::endl;
      }
      //      else{ std::cout << "Didn't match the CA15 jet to an HTT input fatjet with cone size 1.5" << std::endl; }

      theJetCA15HTTTopMass.push_back(theHTTTopMass);
      theJetCA15HTTfrec.push_back(theHTTfrec);

      theNjettinessTau1 = -std::numeric_limits<double>::max();
      theNjettinessTau2 = -std::numeric_limits<double>::max();
      theNjettinessTau3 = -std::numeric_limits<double>::max();
      theNjettinessTau1 = (double)ijet->userFloat("NjettinessCA15CHS:tau1");
      theNjettinessTau2 = (double)ijet->userFloat("NjettinessCA15CHS:tau2");
      theNjettinessTau3 = (double)ijet->userFloat("NjettinessCA15CHS:tau3");
      
      theJetCA15NjettinessTau1.push_back(theNjettinessTau1);
      theJetCA15NjettinessTau2.push_back(theNjettinessTau2);
      theJetCA15NjettinessTau3.push_back(theNjettinessTau3);

    }

    SetValue("theJetCA15Pt", theJetCA15Pt);
    SetValue("theJetCA15Eta", theJetCA15Eta);
    SetValue("theJetCA15Phi", theJetCA15Phi);
    SetValue("theJetCA15Energy", theJetCA15Energy);

    SetValue("theJetCA15HTTTopMass", theJetCA15HTTTopMass);
    SetValue("theJetCA15HTTfrec", theJetCA15HTTfrec);
    
    SetValue("theJetCA15NjettinessTau1", theJetCA15NjettinessTau1);
    SetValue("theJetCA15NjettinessTau2", theJetCA15NjettinessTau2);
    SetValue("theJetCA15NjettinessTau3", theJetCA15NjettinessTau3);
    
    return 0;
}

int JetSubCalc::EndJob()
{
    return 0;
}
