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
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "LJMet/Com/interface/HTTTopJetTagInfo.h"

using namespace std;

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
  edm::InputTag slimmedJetsAK8SDColl_it;
  edm::InputTag slimmedJetsAK8CTTColl_it;
  edm::InputTag selectedPatJetsCA15Coll_it;
  edm::InputTag httTagInfo_it;
  std::string bDiscriminant;
  std::string tagInfo;
  double kappa;
  bool useHTT;
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
    if (mPset.exists("useHTT")) useHTT = mPset.getParameter<bool>("useHTT");
    else useHTT = false;

    if (mPset.exists("slimmedJetColl")){ slimmedJetColl_it = mPset.getParameter<edm::InputTag>("slimmedJetColl");}
    else slimmedJetColl_it = edm::InputTag("slimmedJets");
    
    if (mPset.exists("slimmedJetsAK8Coll")) slimmedJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("slimmedJetsAK8Coll");
    else slimmedJetsAK8Coll_it = edm::InputTag("slimmedJetsAK8");

    if (mPset.exists("bDiscriminant")) bDiscriminant = mPset.getParameter<std::string>("bDiscriminant");
    else bDiscriminant = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
    
    std::cout << " JetSubCalc Bdisc = " << bDiscriminant << std::endl;
    
    if (mPset.exists("tagInfo")) tagInfo = mPset.getParameter<std::string>("tagInfo");
    else tagInfo = "CATop";

    if(mPset.exists("kappa")) kappa = mPset.getParameter<double>("kappa");
    else kappa = 0.5;

    if(useHTT){
      cout << " JetSubCalc is using HTT -- have you installed HEPTopTagger v2?" << endl;

      if (mPset.exists("httTagInfo")) httTagInfo_it = mPset.getParameter<edm::InputTag>("httTagInfo");
      else httTagInfo_it = edm::InputTag("HTT");
  
      if (mPset.exists("selectedJetsCA15Coll")) selectedPatJetsCA15Coll_it = mPset.getParameter<edm::InputTag>("selectedJetsCA15Coll");
      else selectedPatJetsCA15Coll_it = edm::InputTag("selectedPatJetsCA15PFCHSNoHF");
    }
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

    std::vector<double> theJetCEmEnergy;
    std::vector<double> theJetNEmEnergy;
    std::vector<double> theJetCEmEFrac;
    std::vector<double> theJetNEmEFrac;
    std::vector<double> theJetCHadEnergy;
    std::vector<double> theJetNHadEnergy;
    std::vector<double> theJetCHadEFrac;
    std::vector<double> theJetNHadEFrac;
    
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

      theJetCEmEnergy.push_back(theJets[ii]->chargedEmEnergy());
      theJetNEmEnergy.push_back(theJets[ii]->neutralEmEnergy());
      theJetCEmEFrac.push_back(theJets[ii]->chargedEmEnergyFraction());
      theJetNEmEFrac.push_back(theJets[ii]->neutralEmEnergyFraction());	    

      theJetCHadEnergy.push_back(theJets[ii]->chargedHadronEnergy());
      theJetNHadEnergy.push_back(theJets[ii]->neutralHadronEnergy());
      theJetCHadEFrac.push_back(theJets[ii]->chargedHadronEnergyFraction());
      theJetNHadEFrac.push_back(theJets[ii]->neutralHadronEnergyFraction());
      
      theJetIndex.push_back(index);
      theJetnDaughters.push_back((int)theJets[ii]->numberOfDaughters());
 
      //HT
      theJetHT += lv.Pt(); 
    }

    double leading_pt = -999.;
    for (unsigned int ii = 0; ii < theCorrBtagJets.size(); ii++){
        TLorentzVector lv = theCorrBtagJets[ii].first;
        
	if (lv.Pt() > leading_pt){leading_pt = lv.Pt();}
    }
    double second_leading_pt =-999.;
    for (unsigned int ii = 0; ii < theCorrBtagJets.size(); ii++){
        TLorentzVector lv = theCorrBtagJets[ii].first;
        
        if(lv.Pt() > second_leading_pt && lv.Pt() < leading_pt){
	  second_leading_pt = lv.Pt();
        }
    }


    SetValue("theJetPt",     theJetPt);
    SetValue("theJetEta",    theJetEta);
    SetValue("theJetPhi",    theJetPhi);
    SetValue("theJetEnergy", theJetEnergy);
    SetValue("theJetCSV",    theJetCSV);
    SetValue("theJetBTag",   theJetBTag);
    SetValue("theJetFlav",   theJetFlav);

    SetValue("theJetCEmEnergy", theJetCEmEnergy); 
    SetValue("theJetNEmEnergy", theJetNEmEnergy); 
    SetValue("theJetCEmEFrac",  theJetCEmEFrac);  
    SetValue("theJetNEmEFrac",  theJetNEmEFrac);  
    SetValue("theJetCHadEnergy",theJetCHadEnergy);
    SetValue("theJetNHadEnergy",theJetNHadEnergy);
    SetValue("theJetCHadEFrac", theJetCHadEFrac); 
    SetValue("theJetNHadEFrac", theJetNHadEFrac); 

    SetValue("theJetHT", theJetHT);    
    SetValue("theJetLeadPt", leading_pt);
    SetValue("theJetSubLeadPt", second_leading_pt);

    SetValue("theJetVtxMass",     theJetVtxMass);
    SetValue("theJetVtxNtracks",  theJetVtxNtracks);
    SetValue("theJetVtx3DVal",    theJetVtx3DVal);
    SetValue("theJetVtx3DSig",    theJetVtx3DSig);
    SetValue("theJetPileupJetId", theJetPileupJetId);
    SetValue("theJetnDaughters", theJetnDaughters);

    // Load in AK8 jets (no selection performed on these)
    edm::Handle<std::vector<pat::Jet> > theAK8Jets;
    event.getByLabel(slimmedJetsAK8Coll_it, theAK8Jets);

    // Four vector
    std::vector<double> theJetAK8Pt;
    std::vector<double> theJetAK8Eta;
    std::vector<double> theJetAK8Phi;
    std::vector<double> theJetAK8Energy;
    std::vector<double> theJetAK8CSV;
    std::vector<double> theJetAK8JetCharge;

    std::vector<double> theJetAK8CEmEnergy;
    std::vector<double> theJetAK8NEmEnergy;
    std::vector<double> theJetAK8CEmEFrac;
    std::vector<double> theJetAK8NEmEFrac;
    std::vector<double> theJetAK8CHadEnergy;
    std::vector<double> theJetAK8NHadEnergy;
    std::vector<double> theJetAK8CHadEFrac;
    std::vector<double> theJetAK8NHadEFrac;
    
    // Pruned, trimmed and filtered masses available
    std::vector<double> theJetAK8PrunedMass;
    std::vector<double> theJetAK8TrimmedMass;
    std::vector<double> theJetAK8FilteredMass;
    std::vector<double> theJetAK8SoftDropMass;
    
    // n-subjettiness variables tau1, tau2, and tau3 available
    std::vector<double> theJetAK8NjettinessTau1;
    std::vector<double> theJetAK8NjettinessTau2;
    std::vector<double> theJetAK8NjettinessTau3;

    std::vector<bool>   theJetAK8SoftDropTau32Tag;
    std::vector<bool>   theJetAK8SoftDropTau21Tag;
    std::vector<bool>   theJetAK8PrunedTau21Tag;

    std::vector<int>    theJetAK8caTopGoodMatch;    
    std::vector<double> theJetAK8caTopTopMass;
    std::vector<double> theJetAK8caTopMinMass;
    std::vector<int>    theJetAK8caTopnSubJets;
    std::vector<bool>   theJetAK8caTopRun1Tag;

    std::vector<double> theJetAK8Mass;
    std::vector<int>    theJetAK8Index;
    std::vector<int>    theJetAK8nDaughters;
    
    std::vector<double> theJetAK8SDSubjetPt;
    std::vector<double> theJetAK8SDSubjetEta;
    std::vector<double> theJetAK8SDSubjetPhi;
    std::vector<double> theJetAK8SDSubjetMass;
    std::vector<double> theJetAK8SDSubjetCSV;
    std::vector<double> theJetAK8SDSubjetDR;
    std::vector<int> theJetAK8SDSubjetIndex;
    std::vector<int> theJetAK8SDSubjetSize;
    std::vector<int> theJetAK8SDSubjetNCSVL;
    std::vector<int> theJetAK8SDSubjetNCSVM;

    std::vector<double> theJetAK8caTopSubjetPt;
    std::vector<double> theJetAK8caTopSubjetEta;
    std::vector<double> theJetAK8caTopSubjetPhi;
    std::vector<double> theJetAK8caTopSubjetMass;
    std::vector<double> theJetAK8caTopSubjetCSV;
    std::vector<double> theJetAK8caTopSubjetDR;
    std::vector<int> theJetAK8caTopSubjetIndex;
    std::vector<int> theJetAK8caTopSubjetSize;
    std::vector<int> theJetAK8caTopSubjetNCSVL;
    std::vector<int> theJetAK8caTopSubjetNCSVM;
    
    double topMass, minMass, jetCharge;
    int nSubJets;
    double thePrunedMass, theTrimmedMass, theFilteredMass, theSoftDropMass;
    double theNjettinessTau1, theNjettinessTau2, theNjettinessTau3;

    double SDsubjetPt;
    double SDsubjetEta; 
    double SDsubjetPhi;      
    double SDsubjetMass;      
    double SDsubjetBdisc;     
    double SDdeltaRsubjetJet; 
    double CTTsubjetPt;
    double CTTsubjetEta; 
    double CTTsubjetPhi;      
    double CTTsubjetMass;      
    double CTTsubjetBdisc;     
    double CTTdeltaRsubjetJet; 

    int nSDSubJets;
    int nSDSubsCSVL;
    int nSDSubsCSVM;
    int SDSubJetIndex;
    int nTopSubJets;
    int nTopSubsCSVL;
    int nTopSubsCSVM;
    int TopSubJetIndex;

    for (std::vector<pat::Jet>::const_iterator ijet = theAK8Jets->begin(); ijet != theAK8Jets->end(); ijet++) {
      int index = (int)(ijet-theAK8Jets->begin());
      
      thePrunedMass   = -std::numeric_limits<double>::max();
      theTrimmedMass  = -std::numeric_limits<double>::max();
      theFilteredMass = -std::numeric_limits<double>::max();
      theSoftDropMass = -std::numeric_limits<double>::max();
      
      thePrunedMass   = (double)ijet->userFloat("ak8PFJetsCHSPrunedMass");
      theTrimmedMass  = (double)ijet->userFloat("ak8PFJetsCHSTrimmedMass");
      theFilteredMass = (double)ijet->userFloat("ak8PFJetsCHSFilteredMass");
      theSoftDropMass = (double)ijet->userFloat("ak8PFJetsCHSSoftDropMass");
      
      theNjettinessTau1 = -std::numeric_limits<double>::max();
      theNjettinessTau2 = -std::numeric_limits<double>::max();
      theNjettinessTau3 = -std::numeric_limits<double>::max();
      theNjettinessTau1 = (double)ijet->userFloat("NjettinessAK8:tau1");
      theNjettinessTau2 = (double)ijet->userFloat("NjettinessAK8:tau2");
      theNjettinessTau3 = (double)ijet->userFloat("NjettinessAK8:tau3");

      double theNjettinessTau21 = 99;
      double theNjettinessTau32 = 99;
      if(theNjettinessTau1!=2) theNjettinessTau21 = theNjettinessTau2/theNjettinessTau1;
      if(theNjettinessTau2!=2) theNjettinessTau32 = theNjettinessTau3/theNjettinessTau2;

      // Soft Drop + Nsubjettiness tagger
      bool SoftDropTau32Tagged = false;
      if(theSoftDropMass<230 && theSoftDropMass>140 && theNjettinessTau32 <0.65) SoftDropTau32Tagged = true;

      bool SoftDropTau21Tagged = false;
      if(ijet->pt() > 200 && theNjettinessTau21 < 0.65 && theSoftDropMass > 55 && theSoftDropMass < 95){
	SoftDropTau21Tagged = true;
      }
      bool PruningTau21Tagged = false;
      if(ijet->pt() > 200 && theNjettinessTau21 < 0.65 && thePrunedMass > 55 && thePrunedMass < 95){
	PruningTau21Tagged = true;
      }
      
      theJetAK8Pt    .push_back(ijet->pt());
      theJetAK8Eta   .push_back(ijet->eta());
      theJetAK8Phi   .push_back(ijet->phi());
      theJetAK8Energy.push_back(ijet->energy());
      theJetAK8CSV.push_back(ijet->bDiscriminator( bDiscriminant ));

      theJetAK8CEmEnergy.push_back(ijet->chargedEmEnergy());
      theJetAK8NEmEnergy.push_back(ijet->neutralEmEnergy());
      theJetAK8CEmEFrac.push_back(ijet->chargedEmEnergyFraction());
      theJetAK8NEmEFrac.push_back(ijet->neutralEmEnergyFraction());	    

      theJetAK8CHadEnergy.push_back(ijet->chargedHadronEnergy());
      theJetAK8NHadEnergy.push_back(ijet->neutralHadronEnergy());
      theJetAK8CHadEFrac.push_back(ijet->chargedHadronEnergyFraction());
      theJetAK8NHadEFrac.push_back(ijet->neutralHadronEnergyFraction());
      
      theJetAK8Mass.push_back(ijet->mass());
      theJetAK8PrunedMass  .push_back(thePrunedMass);
      theJetAK8TrimmedMass .push_back(theTrimmedMass);
      theJetAK8FilteredMass.push_back(theFilteredMass);
      theJetAK8SoftDropMass.push_back(theSoftDropMass);
      
      theJetAK8NjettinessTau1.push_back(theNjettinessTau1);
      theJetAK8NjettinessTau2.push_back(theNjettinessTau2);
      theJetAK8NjettinessTau3.push_back(theNjettinessTau3);
      
      theJetAK8SoftDropTau32Tag.push_back(SoftDropTau32Tagged);
      theJetAK8SoftDropTau21Tag.push_back(SoftDropTau21Tagged);
      theJetAK8PrunedTau21Tag.push_back(PruningTau21Tagged);

      theJetAK8nDaughters.push_back((int)ijet->numberOfDaughters());
      theJetAK8Index.push_back(index);

      //JetCharge calculation

      reco::Jet::Constituents constituents = ijet->getJetConstituents();

      double sumWeightedCharge = 0.0;
      int con_charge = 0;
      double con_pt = 0.0;
      for(auto constituentItr=constituents.begin(); constituentItr!=constituents.end(); ++constituentItr){
	edm::Ptr<reco::Candidate> constituent=*constituentItr;

	con_charge = (int)constituent->charge();
	con_pt     = (double)constituent->pt();
	
	sumWeightedCharge = sumWeightedCharge + ( con_charge * pow(con_pt,kappa) );
	
      }

      jetCharge  = 1.0/( pow( (ijet->pt()), kappa) ) * sumWeightedCharge;

      theJetAK8JetCharge.push_back(jetCharge);

      // CMSTopTagger information
      bool Run1CMStopTagged = false;
      reco::CATopJetTagInfo const * jetInfo = dynamic_cast<reco::CATopJetTagInfo const *>( ijet->tagInfo( tagInfo ));
      topMass   = -std::numeric_limits<double>::max();
      minMass   = -std::numeric_limits<double>::max();
      nSubJets  = std::numeric_limits<int>::min();

      if(jetInfo != 0){
	topMass = jetInfo->properties().topMass;
	minMass = jetInfo->properties().minMass;
	nSubJets = jetInfo->properties().nSubJets;
	
	if (nSubJets > 2 && minMass > 50.0 && topMass > 140.0 && topMass < 250.0) Run1CMStopTagged = true; 	
	
	theJetAK8caTopTopMass.push_back(topMass);
	theJetAK8caTopMinMass.push_back(minMass);
	theJetAK8caTopnSubJets.push_back(nSubJets);
	theJetAK8caTopRun1Tag.push_back(Run1CMStopTagged);
      }else{
	theJetAK8caTopTopMass.push_back(topMass);
	theJetAK8caTopMinMass.push_back(minMass);
	theJetAK8caTopnSubJets.push_back(nSubJets);
	theJetAK8caTopRun1Tag.push_back(Run1CMStopTagged);
      }

      // Get Soft drop subjets for subjet b-tagging
      SDSubJetIndex = (int)theJetAK8SDSubjetPt.size();
      nSDSubJets  = std::numeric_limits<int>::min();
      nSDSubsCSVL = 0;
      nSDSubsCSVM = 0;

      auto const & sdSubjets = ijet->subjets("SoftDrop");
      nSDSubJets = (int)sdSubjets.size();
      for ( auto const & it : sdSubjets ) {
	SDsubjetPt        = -std::numeric_limits<double>::max();
	SDsubjetEta       = -std::numeric_limits<double>::max();
	SDsubjetPhi       = -std::numeric_limits<double>::max();
	SDsubjetMass      = -std::numeric_limits<double>::max();
	SDsubjetBdisc     = -std::numeric_limits<double>::max();
	SDdeltaRsubjetJet = std::numeric_limits<double>::max();
      
	SDsubjetPt        = it->pt();
	SDsubjetEta       = it->eta();
	SDsubjetPhi       = it->phi();
	SDsubjetMass      = it->mass();
	SDsubjetBdisc     = it->bDiscriminator(bDiscriminant); 
	SDdeltaRsubjetJet = deltaR(ijet->eta(), ijet->phi(), SDsubjetEta, SDsubjetPhi);

	if(SDsubjetBdisc > 0.605) nSDSubsCSVL++;
	if(SDsubjetBdisc > 0.890) nSDSubsCSVM++;

	theJetAK8SDSubjetPt.push_back(SDsubjetPt);
	theJetAK8SDSubjetEta.push_back(SDsubjetEta);
	theJetAK8SDSubjetPhi.push_back(SDsubjetPhi);
	theJetAK8SDSubjetMass.push_back(SDsubjetMass);
	theJetAK8SDSubjetCSV.push_back(SDsubjetBdisc);
	theJetAK8SDSubjetDR.push_back(SDdeltaRsubjetJet);
      }

      theJetAK8SDSubjetIndex.push_back(SDSubJetIndex);
      theJetAK8SDSubjetSize.push_back(nSDSubJets);
      theJetAK8SDSubjetNCSVL.push_back(nSDSubsCSVL);
      theJetAK8SDSubjetNCSVM.push_back(nSDSubsCSVM);

      // Get top tagged subjets 
      TopSubJetIndex = (int)theJetAK8caTopSubjetPt.size();
      nTopSubJets  = std::numeric_limits<int>::min();
      nTopSubsCSVL  = 0;
      nTopSubsCSVM  = 0;

      auto const & topSubjets = ijet->subjets("CMSTopTag");
      nTopSubJets = (int)topSubjets.size();
      for ( auto const & it : topSubjets ) {
	CTTsubjetPt        = -std::numeric_limits<double>::max();
	CTTsubjetEta       = -std::numeric_limits<double>::max();
	CTTsubjetPhi       = -std::numeric_limits<double>::max();
	CTTsubjetMass      = -std::numeric_limits<double>::max();
	CTTsubjetBdisc     = -std::numeric_limits<double>::max();
	CTTdeltaRsubjetJet = std::numeric_limits<double>::max();

	CTTsubjetPt       = it->pt();
	CTTsubjetEta      = it->eta();
	CTTsubjetPhi      = it->phi();
	CTTsubjetMass     = it->mass();
	CTTsubjetBdisc    = it->bDiscriminator(bDiscriminant); 
	CTTdeltaRsubjetJet = deltaR(ijet->eta(), ijet->phi(), CTTsubjetEta, CTTsubjetPhi);

	if(CTTsubjetBdisc > 0.605) nTopSubsCSVL++;
	if(CTTsubjetBdisc > 0.890) nTopSubsCSVM++;

	theJetAK8caTopSubjetPt.push_back(SDsubjetPt);
	theJetAK8caTopSubjetEta.push_back(SDsubjetEta);
	theJetAK8caTopSubjetPhi.push_back(SDsubjetPhi);
	theJetAK8caTopSubjetMass.push_back(SDsubjetMass);
	theJetAK8caTopSubjetCSV.push_back(SDsubjetBdisc);
	theJetAK8caTopSubjetDR.push_back(SDdeltaRsubjetJet);

      }

      theJetAK8caTopSubjetIndex.push_back(TopSubJetIndex);
      theJetAK8caTopSubjetSize.push_back(nTopSubJets);
      theJetAK8caTopSubjetNCSVL.push_back(nTopSubsCSVL);
      theJetAK8caTopSubjetNCSVM.push_back(nTopSubsCSVM);

    }

    SetValue("theJetAK8Pt",     theJetAK8Pt);
    SetValue("theJetAK8Eta",    theJetAK8Eta);
    SetValue("theJetAK8Phi",    theJetAK8Phi);
    SetValue("theJetAK8Energy", theJetAK8Energy);
    SetValue("theJetAK8CSV",    theJetAK8CSV);
    SetValue("theJetAK8JetCharge", theJetAK8JetCharge);

    SetValue("theJetAK8CEmEnergy", theJetAK8CEmEnergy); 
    SetValue("theJetAK8NEmEnergy", theJetAK8NEmEnergy); 
    SetValue("theJetAK8CEmEFrac",  theJetAK8CEmEFrac);  
    SetValue("theJetAK8NEmEFrac",  theJetAK8NEmEFrac);  
    SetValue("theJetAK8CHadEnergy",theJetAK8CHadEnergy);
    SetValue("theJetAK8NHadEnergy",theJetAK8NHadEnergy);
    SetValue("theJetAK8CHadEFrac", theJetAK8CHadEFrac); 
    SetValue("theJetAK8NHadEFrac", theJetAK8NHadEFrac); 

    SetValue("theJetAK8PrunedMass",   theJetAK8PrunedMass);
    SetValue("theJetAK8TrimmedMass",  theJetAK8TrimmedMass);
    SetValue("theJetAK8FilteredMass", theJetAK8FilteredMass);
    SetValue("theJetAK8SoftDropMass", theJetAK8SoftDropMass);
    
    SetValue("theJetAK8NjettinessTau1", theJetAK8NjettinessTau1);
    SetValue("theJetAK8NjettinessTau2", theJetAK8NjettinessTau2);
    SetValue("theJetAK8NjettinessTau3", theJetAK8NjettinessTau3);

    SetValue("theJetAK8SoftDropTau32Tag", theJetAK8SoftDropTau32Tag);
    SetValue("theJetAK8SoftDropTau21Tag", theJetAK8SoftDropTau21Tag);
    SetValue("theJetAK8PrunedTau21Tag", theJetAK8PrunedTau21Tag);
    
    SetValue("theJetAK8Mass",   theJetAK8Mass);
    SetValue("theJetAK8nDaughters", theJetAK8nDaughters);
    
    SetValue("theJetAK8caTopTopMass", theJetAK8caTopTopMass);
    SetValue("theJetAK8caTopMinMass", theJetAK8caTopMinMass);
    SetValue("theJetAK8caTopnSubJets", theJetAK8caTopnSubJets);
    SetValue("theJetAK8caTopRun1Tag", theJetAK8caTopRun1Tag);

    SetValue("theJetAK8SDSubjetPt",   theJetAK8SDSubjetPt);   
    SetValue("theJetAK8SDSubjetEta",  theJetAK8SDSubjetEta);  
    SetValue("theJetAK8SDSubjetPhi",  theJetAK8SDSubjetPhi);  
    SetValue("theJetAK8SDSubjetMass", theJetAK8SDSubjetMass); 
    SetValue("theJetAK8SDSubjetCSV",  theJetAK8SDSubjetCSV);  
    SetValue("theJetAK8SDSubjetDR",   theJetAK8SDSubjetDR);   
    SetValue("theJetAK8SDSubjetIndex",theJetAK8SDSubjetIndex);
    SetValue("theJetAK8SDSubjetSize", theJetAK8SDSubjetSize); 
    SetValue("theJetAK8SDSubjetNCSVL",theJetAK8SDSubjetNCSVL);
    SetValue("theJetAK8SDSubjetNCSVM",theJetAK8SDSubjetNCSVM);

    SetValue("theJetAK8caTopSubjetPt",   theJetAK8caTopSubjetPt);   
    SetValue("theJetAK8caTopSubjetEta",	 theJetAK8caTopSubjetEta);	 
    SetValue("theJetAK8caTopSubjetPhi",	 theJetAK8caTopSubjetPhi);	 
    SetValue("theJetAK8caTopSubjetMass", theJetAK8caTopSubjetMass); 
    SetValue("theJetAK8caTopSubjetCSV",	 theJetAK8caTopSubjetCSV);	 
    SetValue("theJetAK8caTopSubjetDR",	 theJetAK8caTopSubjetDR);	 
    SetValue("theJetAK8caTopSubjetIndex",theJetAK8caTopSubjetIndex);
    SetValue("theJetAK8caTopSubjetSize", theJetAK8caTopSubjetSize); 
    SetValue("theJetAK8caTopSubjetNCSVL",theJetAK8caTopSubjetNCSVL);
    SetValue("theJetAK8caTopSubjetNCSVM",theJetAK8caTopSubjetNCSVM);

    if(useHTT){

      ///////////////////////////////////////////////////////////////////////
      ///  
      ///   NOTE: You must install HTT to load the HTTTopJetTagInfo object!
      ///   In your CMSSW_X_X_X/src/ directory:
      ///   git cms-merge-topic gkasieczka:htt-v2-74X
      ///   scramv1 b -r -j 8
      ///   
      ///////////////////////////////////////////////////////////////////////
    
      std::vector<double> theJetCA15Pt;
      std::vector<double> theJetCA15Eta;
      std::vector<double> theJetCA15Phi;
      std::vector<double> theJetCA15Energy;
      std::vector<double> theJetCA15Mass;
      std::vector<double> theJetCA15SoftDropMass;
      std::vector<double> theJetCA15NjettinessTau1;
      std::vector<double> theJetCA15NjettinessTau2;
      std::vector<double> theJetCA15NjettinessTau3;
      
      std::vector<double> theJetHTTTopMass;
      std::vector<double> theJetHTTfRec;
      std::vector<double> theJetHTTRopt;
      std::vector<double> theJetHTTRoptCalc;
      std::vector<double> theJetHTTRoptCalcPt;
      std::vector<double> theJetHTTNjettinessTau1;
      std::vector<double> theJetHTTNjettinessTau2;
      std::vector<double> theJetHTTNjettinessTau3;
      std::vector<double> theJetHTTufNjettinessTau1;
      std::vector<double> theJetHTTufNjettinessTau2;
      std::vector<double> theJetHTTufNjettinessTau3;

      double theHTTTopMass;
      double theHTTfRec;
      double theHTTRopt;
      double theHTTRoptCalc;
      double theHTTRoptCalcPt;
      double NjettinessTau1;
      double NjettinessTau2;
      double NjettinessTau3;
      
      edm::Handle<std::vector<pat::Jet> > theCA15Jets;
      event.getByLabel(selectedPatJetsCA15Coll_it, theCA15Jets);

      cout << "Created CA15 Jets, getting HTT" << endl;
      edm::Handle<std::vector<reco::HTTTopJetTagInfo> > httTagInfo;
      cout << "Created HTT handle" << endl;
      event.getByLabel(httTagInfo_it, httTagInfo);
      cout << "Got HTT collection" << endl;
      for (std::vector<pat::Jet>::const_iterator ijet = theCA15Jets->begin(); ijet != theCA15Jets->end(); ijet++){

	theJetCA15Pt.push_back(ijet->pt());
	theJetCA15Eta.push_back(ijet->eta());
	theJetCA15Phi.push_back(ijet->phi());
	theJetCA15Energy.push_back(ijet->energy());
	theJetCA15Mass.push_back(ijet->mass());

	theNjettinessTau1 = -std::numeric_limits<double>::max();
	theNjettinessTau2 = -std::numeric_limits<double>::max();
	theNjettinessTau3 = -std::numeric_limits<double>::max();
	theNjettinessTau1 = (double)ijet->userFloat("NjettinessCA15CHSNoHF:tau1");
	theNjettinessTau2 = (double)ijet->userFloat("NjettinessCA15CHSNoHF:tau2");
	theNjettinessTau3 = (double)ijet->userFloat("NjettinessCA15CHSNoHF:tau3");
      
	theJetCA15NjettinessTau1.push_back(theNjettinessTau1);
	theJetCA15NjettinessTau2.push_back(theNjettinessTau2);
	theJetCA15NjettinessTau3.push_back(theNjettinessTau3);
	
	theSoftDropMass = -std::numeric_limits<double>::max();
	theSoftDropMass = (double)ijet->userFloat("ca15PFJetsCHSNoHFSoftDropMass");
	theJetCA15SoftDropMass.push_back(theSoftDropMass);
	
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
	
 	theHTTTopMass = -99.99;
	theHTTfRec = -99.99;
	theHTTRopt = -99.99;
	theHTTRoptCalc = -99.99;
	theHTTRoptCalcPt = -99.99;
	theNjettinessTau1 = -99.99;
	theNjettinessTau2 = -99.99;
	theNjettinessTau3 = -99.99;
	NjettinessTau1 = -99.99;
	NjettinessTau2 = -99.99;
	NjettinessTau3 = -99.99;

	if(minDR < 1.5){
	  theHTTTopMass = CA15HTTProperties.topMass;    // mass of the top quark candidate
	  theHTTfRec = CA15HTTProperties.fRec;          // min distance of m_ij/m_123 from m_W/m_top
	  theHTTRopt = CA15HTTProperties.Ropt;          // R_opt found in Optimal procedure
	  theHTTRoptCalc = CA15HTTProperties.RoptCalc;     // R_opt calc for a top based on filtered fat jet pT
	  theHTTRoptCalcPt = CA15HTTProperties.ptForRoptCalc; // Filtered fat jet pT for R_opt calculation
	  theNjettinessTau1 = CA15HTTProperties.tau1Filtered; 
	  theNjettinessTau2 = CA15HTTProperties.tau2Filtered;
	  theNjettinessTau3 = CA15HTTProperties.tau3Filtered;
	  NjettinessTau1 = CA15HTTProperties.tau1Unfiltered;
	  NjettinessTau2 = CA15HTTProperties.tau2Unfiltered;
	  NjettinessTau3 = CA15HTTProperties.tau3Unfiltered;
	}
	
	theJetHTTTopMass.push_back(theHTTTopMass);
	theJetHTTfRec.push_back(theHTTfRec);
	theJetHTTRopt.push_back(theHTTRopt);
	theJetHTTRoptCalc.push_back(theHTTRoptCalc);
	theJetHTTRoptCalcPt.push_back(theHTTRoptCalcPt);
	theJetHTTNjettinessTau1.push_back(theNjettinessTau1);
	theJetHTTNjettinessTau2.push_back(theNjettinessTau2);
	theJetHTTNjettinessTau3.push_back(theNjettinessTau3);
	theJetHTTufNjettinessTau1.push_back(NjettinessTau1);
	theJetHTTufNjettinessTau2.push_back(NjettinessTau2);
	theJetHTTufNjettinessTau3.push_back(NjettinessTau3);
	
      }
      
      SetValue("theJetCA15Pt", theJetCA15Pt);
      SetValue("theJetCA15Eta", theJetCA15Eta);
      SetValue("theJetCA15Phi", theJetCA15Phi);
      SetValue("theJetCA15Energy", theJetCA15Energy);
      SetValue("theJetCA15Mass", theJetCA15Mass);
      SetValue("theJetCA15SoftDropMass", theJetCA15SoftDropMass);
      SetValue("theJetCA15NjettinessTau1", theJetCA15NjettinessTau1);
      SetValue("theJetCA15NjettinessTau2", theJetCA15NjettinessTau2);
      SetValue("theJetCA15NjettinessTau3", theJetCA15NjettinessTau3);
      
      SetValue("theJetHTTTopMass", theJetHTTTopMass);
      SetValue("theJetHTTfRec", theJetHTTfRec);
      SetValue("theJetHTTRopt", theJetHTTRopt);
      SetValue("theJetHTTRoptCalc", theJetHTTRoptCalc);
      SetValue("theJetHTTRoptCalcPt", theJetHTTRoptCalcPt);
      
      SetValue("theJetHTTNjettinessTau1",theJetHTTNjettinessTau1);    
      SetValue("theJetHTTNjettinessTau2",theJetHTTNjettinessTau2);
      SetValue("theJetHTTNjettinessTau3",theJetHTTNjettinessTau3);
      SetValue("theJetHTTufNjettinessTau1",theJetHTTufNjettinessTau1);    
      SetValue("theJetHTTufNjettinessTau2",theJetHTTufNjettinessTau2);
      SetValue("theJetHTTufNjettinessTau3",theJetHTTufNjettinessTau3);
    }

    return 0;
}

int JetSubCalc::EndJob()
{
    return 0;
}
