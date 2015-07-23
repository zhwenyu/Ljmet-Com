/*
  Calculator for the Run 2 inclusive TprimeTprime analysis
  Finds Tprime particles and tags the decays chain so branching ratios can be scanned

  Author: Julie Hogan, 2015
*/


#include <iostream>
#include <algorithm>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h" 
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h"
class LjmetFactory;


class TpTpCalc : public BaseCalc{
  
public:
  
    TpTpCalc();
    virtual ~TpTpCalc(){}

    virtual int BeginJob(){
 
      if (mPset.exists("genParticles")) genParticles_it = mPset.getParameter<edm::InputTag>("genParticles");
      else                              genParticles_it = edm::InputTag("prunedGenParticles");
 
      return 0;
    }

    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}

  
private:
  
  edm::InputTag genParticles_it;
 
};



static int reg = LjmetFactory::GetInstance()->Register(new TpTpCalc(), "TpTpCalc");



TpTpCalc::TpTpCalc(){
}

int TpTpCalc::AnalyzeEvent(edm::EventBase const & event,
                             BaseEventSelector * selector){
    //
    // compute event variables here
    //
      
    bool isTZTZ = false;
    bool isTZTH = false;
    bool isTZBW = false;
    bool isTHTH = false;
    bool isTHBW = false;
    bool isBWBW = false;

    std::vector<int>    tPrimeStatus;
    std::vector<int>    tPrimeID;
    std::vector<int>    tPrimeNDaughters;
    std::vector<double> tPrimeMass;
    std::vector<double> tPrimePt;
    std::vector<double> tPrimeEta;
    std::vector<double> tPrimePhi;
    std::vector<double> tPrimeEnergy;

    std::vector<int>    quarkID;
    std::vector<double> quarkPt;
    std::vector<double> quarkEta;
    std::vector<double> quarkPhi;
    std::vector<double> quarkEnergy;

    std::vector<int>    bosonID;
    std::vector<double> bosonPt;
    std::vector<double> bosonEta;
    std::vector<double> bosonPhi;
    std::vector<double> bosonEnergy;

    // Get the generated particle collection
    edm::Handle<reco::GenParticleCollection> genParticles;
    event.getByLabel(genParticles_it, genParticles);

    // loop over all gen particles in event
    for(size_t i = 0; i < genParticles->size(); i++){
      const reco::GenParticle &p = (*genParticles).at(i);
      int id = p.pdgId();

      // find Tprime particles (+/- 8)
      if(abs(id) != 8) continue;

      tPrimeStatus.push_back( p.status() );
      tPrimeID.push_back( p.pdgId() );
      tPrimeMass.push_back( p.mass() );
      tPrimePt.push_back( p.pt() );
      tPrimeEta.push_back( p.eta() );
      tPrimePhi.push_back( p.phi() );
      tPrimeEnergy.push_back( p.energy() );
      tPrimeNDaughters.push_back( p.numberOfDaughters() );
      
      // get the Tprime daughter particles (looking for t, b, Z, H, W)
      size_t nDs = p.numberOfDaughters();
      for(size_t i = 0; i < nDs; i++){
	int dauId = (p.daughter(i))->pdgId();
	const reco::Candidate *d = p.daughter(i);
	if(d->pdgId() != dauId) std::cout << "making daughter GenParticle didn't work" << std::endl;
	
	// start with the tops and bottoms
	if(abs(dauId) == 5 || abs(dauId) == 6){
	  quarkID.push_back(dauId);
	  quarkPt.push_back(d->pt());
	  quarkEta.push_back(d->eta());
	  quarkPhi.push_back(d->phi());
	  quarkEnergy.push_back(d->energy());
	}
	
	// next are the electroweak bosons
	else if(abs(dauId) > 22 && abs(dauId) < 26){
	  bosonID.push_back(dauId);
	  bosonPt.push_back(d->pt());
	  bosonEta.push_back(d->eta());
	  bosonPhi.push_back(d->phi());
	  bosonEnergy.push_back(d->energy());
	}
	else continue;
      }
    }

    // store variables into the tree
    SetValue("tPrimeStatus",tPrimeStatus);
    SetValue("tPrimeID",tPrimeID);
    SetValue("tPrimeMass",tPrimeMass);
    SetValue("tPrimePt",tPrimePt);
    SetValue("tPrimeEta",tPrimeEta);
    SetValue("tPrimePhi",tPrimePhi);
    SetValue("tPrimeEnergy",tPrimeEnergy);
    SetValue("tPrimeNDaughters",tPrimeNDaughters);

    SetValue("quarkID",quarkID);
    SetValue("quarkPt",quarkPt);
    SetValue("quarkEta",quarkEta);
    SetValue("quarkPhi",quarkPhi);
    SetValue("quarkEnergy",quarkEnergy);

    SetValue("bosonID",bosonID);
    SetValue("bosonPt",bosonPt);
    SetValue("bosonEta",bosonEta);
    SetValue("bosonPhi",bosonPhi);
    SetValue("bosonEnergy",bosonEnergy);

    // Now time for some checking -- do we have the right particles? 
    if(quarkID.size() == 0){ return 0;}
    else if(quarkID.size() != 2){
      std::cout << "More/less than 2 quarks stored: " << quarkID.size() << std::endl;
      for(size_t i = 0; i < quarkID.size(); i++){std::cout << "quark " << i << " = " << quarkID[i] << std::endl;}

      // First check for opposite sign particles in the first two slots
      // A b-bbar pair can sneak in between the daughters we want
      double test = quarkID[0]*quarkID[1];
      int sign = -1; if(test > 0) sign = 1;

      if(sign > 0){
	if(quarkID.size() == 4){
	  std::swap(quarkID[2],quarkID[3]);
	  std::swap(quarkPt[2],quarkPt[3]);
	  std::swap(quarkEta[2],quarkEta[3]);
	  std::swap(quarkPhi[2],quarkPhi[3]);
	  std::swap(quarkEnergy[2],quarkEnergy[3]);

	  std::cout << "Moved up an opposite sign quark: 3 --> 2" << std::endl;
	}
	std::swap(quarkID[1],quarkID[2]);
	std::swap(quarkPt[1],quarkPt[2]);
	std::swap(quarkEta[1],quarkEta[2]);
	std::swap(quarkPhi[1],quarkPhi[2]);
	std::swap(quarkEnergy[1],quarkEnergy[2]);
	
	std::cout << "Moved up an opposite sign quark: 2 --> 1" << std::endl;

	test = quarkID[0]*quarkID[1];
	sign = -1; if(test > 0) sign = 1;
	if(sign < 0) std::cout << "Signs are fixed!" << std::endl;
      }

      // Now check for t's that come after a b-bbar pair, even if signs are ok
      if(quarkID.size() > 3 && abs(quarkID[3]) == 6){
	std::swap(quarkID[2],quarkID[3]);
	std::swap(quarkPt[2],quarkPt[3]);
	std::swap(quarkEta[2],quarkEta[3]);
	std::swap(quarkPhi[2],quarkPhi[3]);
	std::swap(quarkEnergy[2],quarkEnergy[3]);

	std::cout << "Moved up a top, 3 --> 2" << std::endl;
      }
      if(quarkID.size() > 2 && abs(quarkID[2]) == 6){
	std::swap(quarkID[1],quarkID[2]);
	std::swap(quarkPt[1],quarkPt[2]);
	std::swap(quarkEta[1],quarkEta[2]);
	std::swap(quarkPhi[1],quarkPhi[2]);
	std::swap(quarkEnergy[1],quarkEnergy[2]);

	std::cout << "Moved up a top, 2 --> 1" << std::endl;
      }
    }
    else{}

    if(bosonID.size() == 0){return 0;}
    else if(bosonID.size() != 2) std::cout << "More/less than 2 bosons stored: " << bosonID.size() << std::endl;
    else{}

    // tag the decay chains according to ID'd quarks and bosons.
    // After the fixes above there should not be errors (if we've done it right!)

    // 2 b quarks, check for matching W's
    if(abs(quarkID[0]) == 5 && abs(quarkID[1]) == 5){
      if(abs(bosonID[0]) == 24 && abs(bosonID[1]) == 24) isBWBW = true;
      else std::cout << "2 b daughters didn't match bWbW: " << bosonID[0] << ", " << bosonID[1] << std::endl;
    }
    // 2 t quarks, check for Z's and H's
    else if(abs(quarkID[0]) == 6 && abs(quarkID[1]) == 6){
      if(bosonID[0] == 23 && bosonID[1] == 23) isTZTZ = true;
      else if(bosonID[0] == 25 && bosonID[1] == 25) isTHTH = true;
      else if(bosonID[0] == 25 && bosonID[1] == 23) isTZTH = true;
      else if(bosonID[0] == 23 && bosonID[1] == 25) isTZTH = true;
      else std::cout << "2 t daughters didn't match tZtZ, tHtH, or tZtH" << bosonID[0] << ", " << bosonID[1] << std::endl;
    }
    // t-b pairs, check for correlating bosons in the right spots
    else if(abs(quarkID[0]) == 6 && abs(quarkID[1]) == 5){
      if(bosonID[0] == 23 && abs(bosonID[1]) == 24) isTZBW = true;
      else if(bosonID[0] == 25 && abs(bosonID[1]) == 24) isTHBW = true;
      else std::cout << "t - b pair didn't match Z/H - W pair" << bosonID[0] << ", " << bosonID[1] << std::endl; 
    }
    // b-t pairs, check for correlating bosons in the right spots
    else if(abs(quarkID[1]) == 6 && abs(quarkID[0]) == 5){
      if(bosonID[1] == 23 && abs(bosonID[0] == 24)) isTZBW = true;
      else if(bosonID[1] == 25 && abs(bosonID[0]) == 24) isTHBW = true;
      else std::cout << "b - t pair didn't match W - Z/H pair" << bosonID[0] << ", " << bosonID[1] << std::endl;
    }
    // error messages if we found something else entirely
    else{
      std::cout << "daughters didn't match a recognized pattern" << std::endl;
      for(size_t i = 0; i < quarkID.size(); i++){ std::cout << "quark " << i << " = " << quarkID[i] << std::endl; }
      for(size_t i = 0; i < bosonID.size(); i++){ std::cout << "boson " << i << " = " << bosonID[i] << std::endl; }
    }
    
    // store tag values in the ROOT tree
    SetValue("isTZTZ",isTZTZ);
    SetValue("isTZTH",isTZTH);
    SetValue("isTZBW",isTZBW);
    SetValue("isTHTH",isTHTH);
    SetValue("isTHBW",isTHBW);
    SetValue("isBWBW",isBWBW);
    /*
    vector<double> genMuPt;
    vector<double> genMuPhi;
    vector<double> genMuEta;
    vector<bool> genMuPrompt;
    vector<bool> genMuHadron;

    // loop over all gen particles in event
    for(size_t i = 0; i < genParticles->size(); i++){
      const reco::GenParticle &p = (*genParticles).at(i);
      int id = p.pdgId();

      // find Tprime particles (+/- 8)
      if(abs(id) != 13) continue;

      genMuPt.push_back(p.pt());
      genMuPhi.push_back(p.phi());
      genMuPhi.push_back(p.eta());
      genMuPrompt.push_back(p.isPrompt());
      genMuHadron.push_back(p.isDirectHadronDecayProduct());
    }

    SetValue("genMuPt",genMuPt);
    SetValue("genMuPhi",genMuPhi);
    SetValue("genMuEta",genMuEta);
    SetValue("genMuHadron",genMuHadron);
    SetValue("genMuPrompt",genMuPrompt);
    */	  
    return 0;
}

