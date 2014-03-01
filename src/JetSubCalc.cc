/*
  Calculator for substructure variables
  
  Author: Joshua Swanson, 2014
*/



#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/Nsubjettiness.hh"
#include "LJMet/Com/interface/Njettiness.hh"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "DataFormats/PatCandidates/interface/Jet.h"


#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"

using std::cout;
using std::endl;

class LjmetFactory;

class JetSubCalc : public BaseCalc{
  
public:
  
  JetSubCalc();
  virtual ~JetSubCalc(){}
  
  virtual int BeginJob();
  virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
  virtual int EndJob(){return 0;}

private:
    

};


static int reg = LjmetFactory::GetInstance()->Register(new JetSubCalc(), "JetSubCalc");


JetSubCalc::JetSubCalc(){
}

int JetSubCalc::BeginJob(){
  return 0;
}

int JetSubCalc::AnalyzeEvent(edm::EventBase const & event,
			       BaseEventSelector * selector){

  
    //Get Top-like jets
    edm::InputTag topJetColl = edm::InputTag("goodPatJetsCATopTagPFPacked");
    edm::Handle<std::vector<pat::Jet> > topJets;
    event.getByLabel(topJetColl, topJets);

    //Four vector
    std::vector <double> CATopJetPt;
    std::vector <double> CATopJetEta;
    std::vector <double> CATopJetPhi;
    std::vector <double> CATopJetEnergy;

    std::vector <double> CATopJetCSV;
    //std::vector <double> CATopJetRCN;

    //Identity
    std::vector <int> CATopJetIndex;
    std::vector <int> CATopJetnDaughters;

    //Top-like properties
    std::vector <double> CATopJetTopMass;
    std::vector <double> CATopJetMinPairMass;

    //Daughter four vector and index
    std::vector <double> CATopDaughterPt;
    std::vector <double> CATopDaughterEta;
    std::vector <double> CATopDaughterPhi;
    std::vector <double> CATopDaughterEnergy;

    std::vector <int> CATopDaughterMotherIndex;

    for (std::vector<pat::Jet>::const_iterator ijet = topJets->begin(); ijet != topJets->end(); ijet++){

      int index = (int)(ijet-topJets->begin());

      CATopJetPt     . push_back(ijet->pt());
      CATopJetEta    . push_back(ijet->eta());
      CATopJetPhi    . push_back(ijet->phi());
      CATopJetEnergy . push_back(ijet->energy());

      CATopJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
      //CATopJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    

      CATopJetIndex      . push_back(index);
      CATopJetnDaughters . push_back((int)ijet->numberOfDaughters());

      reco::CATopJetTagInfo* jetInfo = (reco::CATopJetTagInfo*) ijet->tagInfo("CATop");

      CATopJetTopMass     . push_back(jetInfo->properties().topMass);
      CATopJetMinPairMass . push_back(jetInfo->properties().minMass);

      for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
		CATopDaughterPt     . push_back(ijet->daughter(ui)->pt());
		CATopDaughterEta    . push_back(ijet->daughter(ui)->eta());
		CATopDaughterPhi    . push_back(ijet->daughter(ui)->phi());
		CATopDaughterEnergy . push_back(ijet->daughter(ui)->energy());        

		CATopDaughterMotherIndex . push_back(index);      
      }
    }

    //Four vector
    SetValue("CATopJetPt"     , CATopJetPt);
    SetValue("CATopJetEta"    , CATopJetEta);
    SetValue("CATopJetPhi"    , CATopJetPhi);
    SetValue("CATopJetEnergy" , CATopJetEnergy);

    SetValue("CATopJetCSV"    , CATopJetCSV);
    //SetValue("CATopJetRCN"    , CATopJetRCN);

    //Identity
    SetValue("CATopJetIndex"      , CATopJetIndex);
    SetValue("CATopJetnDaughters" , CATopJetnDaughters);

    //Properties
    SetValue("CATopJetTopMass"     , CATopJetTopMass);
    SetValue("CATopJetMinPairMass" , CATopJetMinPairMass);

    //Daughter four vector and index
    SetValue("CATopDaughterPt"     , CATopDaughterPt);
    SetValue("CATopDaughterEta"    , CATopDaughterEta);
    SetValue("CATopDaughterPhi"    , CATopDaughterPhi);
    SetValue("CATopDaughterEnergy" , CATopDaughterEnergy);

    SetValue("CATopDaughterMotherIndex"      , CATopDaughterMotherIndex);

    //Get CA8 jets for W's
    edm::InputTag CAWJetColl = edm::InputTag("goodPatJetsCA8PrunedPFPacked");
    edm::Handle<std::vector<pat::Jet> > CAWJets;
    event.getByLabel(CAWJetColl, CAWJets);

    //Four vector
    std::vector <double> CAWJetPt;
    std::vector <double> CAWJetEta;
    std::vector <double> CAWJetPhi;
    std::vector <double> CAWJetEnergy;

    std::vector <double> CAWJetCSV;
    //std::vector <double> CAWJetRCN;

    //Identity
    std::vector <int> CAWJetIndex;
    std::vector <int> CAWJetnDaughters;

    //Mass
    std::vector <double> CAWJetMass;

    //Daughter four vector and index
    std::vector <double> CAWDaughterPt;
    std::vector <double> CAWDaughterEta;
    std::vector <double> CAWDaughterPhi;
    std::vector <double> CAWDaughterEnergy;

    std::vector <int> CAWDaughterMotherIndex;

    for (std::vector<pat::Jet>::const_iterator ijet = CAWJets->begin(); ijet != CAWJets->end(); ijet++){

      int index = (int)(ijet-CAWJets->begin());

      //Four vector
      CAWJetPt     . push_back(ijet->pt());
      CAWJetEta    . push_back(ijet->eta());
      CAWJetPhi    . push_back(ijet->phi());
      CAWJetEnergy . push_back(ijet->energy());        

      CAWJetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
      //CAWJetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));    

      //Identity
      CAWJetIndex      . push_back(index);
      CAWJetnDaughters . push_back((int)ijet->numberOfDaughters());

      //Mass
      CAWJetMass . push_back(ijet->mass());

      for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
		CAWDaughterPt     . push_back(ijet->daughter(ui)->pt());
		CAWDaughterEta    . push_back(ijet->daughter(ui)->eta());
		CAWDaughterPhi    . push_back(ijet->daughter(ui)->phi());
		CAWDaughterEnergy . push_back(ijet->daughter(ui)->energy());        

		CAWDaughterMotherIndex . push_back(index);      
      }
	}

    //Four vector
    SetValue("CAWJetPt"     , CAWJetPt);
    SetValue("CAWJetEta"    , CAWJetEta);
    SetValue("CAWJetPhi"    , CAWJetPhi);
    SetValue("CAWJetEnergy" , CAWJetEnergy);

    SetValue("CAWJetCSV"    , CAWJetCSV);
    //SetValue("CAWJetRCN"    , CAWJetRCN);

    //Identity
    SetValue("CAWJetIndex"      , CAWJetIndex);
    SetValue("CAWJetnDaughters" , CAWJetnDaughters);

    //Mass
    SetValue("CAWJetMass"     , CAWJetMass);

    //Daughter four vector and index
    SetValue("CAWDaughterPt"     , CAWDaughterPt);
    SetValue("CAWDaughterEta"    , CAWDaughterEta);
    SetValue("CAWDaughterPhi"    , CAWDaughterPhi);
    SetValue("CAWDaughterEnergy" , CAWDaughterEnergy);

    SetValue("CAWDaughterMotherIndex" , CAWDaughterMotherIndex);

    //Get all CA8 jets (not just for W and Top)
    edm::InputTag CA8JetColl = edm::InputTag("goodPatJetsCA8PF");
    edm::Handle<std::vector<pat::Jet> > CA8Jets;
    event.getByLabel(CA8JetColl, CA8Jets);

    //Four vector
    std::vector <double> CA8JetPt;
    std::vector <double> CA8JetEta;
    std::vector <double> CA8JetPhi;
    std::vector <double> CA8JetEnergy;

    std::vector <double> CA8JetCSV;
    //std::vector <double> CA8JetRCN;
    std::vector <double> CA8Tau1;
    std::vector <double> CA8Tau2;
    std::vector <double> CA8Tau3;
    std::vector <double> CA8Tau4;
    

  	// ---------------------------------------------------------------------------------------------------------
  	// Setup Nsubjettiness
  	// ---------------------------------------------------------------------------------------------------------

  	Nsubjettiness Nsubonepass1(1, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);
  	Nsubjettiness Nsubonepass2(2, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);
  	Nsubjettiness Nsubonepass3(3, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);
  	Nsubjettiness Nsubonepass4(4, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);

    for (std::vector<pat::Jet>::const_iterator ijet = CA8Jets->begin(); ijet != CA8Jets->end(); ijet++){

      //Four vector
      CA8JetPt     . push_back(ijet->pt());
      CA8JetEta    . push_back(ijet->eta());
      CA8JetPhi    . push_back(ijet->phi());
      CA8JetEnergy . push_back(ijet->energy());

      CA8JetCSV    . push_back(ijet->bDiscriminator( "combinedSecondaryVertexBJetTags"));
      //CA8JetRCN    . push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));
      std::vector<fastjet::PseudoJet> FJparticles;
      for (unsigned i = 0; i < ijet->numberOfDaughters() ; i++){
      	const reco::PFCandidate* this_constituent = dynamic_cast<const reco::PFCandidate*>(ijet->daughter(i));
      	FJparticles.push_back( fastjet::PseudoJet( this_constituent->px(),
             	this_constituent->py(),
             	this_constituent->pz(),
             	this_constituent->energy() ) );
      }
	  
	  fastjet::PseudoJet combJet = fastjet::join(FJparticles);

	  CA8Tau1.push_back( Nsubonepass1.result(combJet) );
      CA8Tau2.push_back( Nsubonepass2.result(combJet) );
      CA8Tau3.push_back( Nsubonepass3.result(combJet) );
      CA8Tau4.push_back( Nsubonepass4.result(combJet) );
    
    }

    //Four vector
    SetValue("CA8JetPt"     , CA8JetPt);
    SetValue("CA8JetEta"    , CA8JetEta);
    SetValue("CA8JetPhi"    , CA8JetPhi);
    SetValue("CA8JetEnergy" , CA8JetEnergy);

    SetValue("CA8JetCSV"    , CA8JetCSV);
    //SetValue("CA8JetRCN"    , CA8JetRCN);
    SetValue("CA8Tau1" , CA8Tau1);
    SetValue("CA8Tau2" , CA8Tau2);
    SetValue("CA8Tau3" , CA8Tau3);
    SetValue("CA8Tau4" , CA8Tau4);

  	return 0;

}