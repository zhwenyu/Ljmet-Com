/*
  Calculator for the Run 2 inclusive TprimeTprime analysis
  Finds Tprime or Bprime particles and tags the decays chain so branching ratios can be scanned

  Author: Julie Hogan, 2015
*/


#include <iostream>
#include <algorithm>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h" 
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

class LjmetFactory;


class TTbarMassCalc : public BaseCalc{
  
public:
  
    TTbarMassCalc();
    virtual ~TTbarMassCalc(){}

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



//static int reg = LjmetFactory::GetInstance()->Register(new TTbarMassCalc(), "TTbarMassCalc");



TTbarMassCalc::TTbarMassCalc(){
}

int TTbarMassCalc::AnalyzeEvent(edm::EventBase const & event,
                             BaseEventSelector * selector){
    //
    // compute event variables here
    //
      
    std::vector<int>    topID;
    std::vector<int>    topMotherID;
    std::vector<double> topPt;
    std::vector<double> topEta;
    std::vector<double> topPhi;
    std::vector<double> topEnergy;
    std::vector<double> topMass;

    std::vector<int>    allTopsID;
    std::vector<int>    allTopsStatus;
    std::vector<double> allTopsPt;
    std::vector<double> allTopsEta;
    std::vector<double> allTopsPhi;
    std::vector<double> allTopsEnergy;

    double ttbarMass = -999.9;
    TLorentzVector top1;
    TLorentzVector top2;

    // Get the generated particle collection
    edm::Handle<reco::GenParticleCollection> genParticles;
    if(event.getByLabel(genParticles_it, genParticles)){

      // loop over all gen particles in event
      for(size_t i = 0; i < genParticles->size(); i++){
	const reco::GenParticle &p = (*genParticles).at(i);
	int id = p.pdgId();
	
	// find tops 
	if(abs(id) != 6) continue;
	if(abs(p.mother()->pdgId()) == 6) continue;      
	
	topID.push_back(id);
	topPt.push_back(p.pt());
	topEta.push_back(p.eta());
	topPhi.push_back(p.phi());
	topMass.push_back(p.mass());
	topEnergy.push_back(p.energy());
	topMotherID.push_back(p.mother()->pdgId());

	//	cout << "status = " << p.status() << endl;
	
      }
            
      // Now time for some checking -- do we have the right particles? 
      if(topID.size() == 0) ttbarMass = -999;
      else if(topID.size() != 2){
	//	std::cout << "TTbarMassCalc: More/less than 2 quarks stored: " << topID.size() << std::endl;
	//	for(size_t i = 0; i < topID.size(); i++){	
	//	  std::cout << "quark " << i << " = " << topID[i] << ", mother ID = " << topMotherID[i] << std::endl;
	//	}
      }
      else{
	if(topID[0]*topID[1] > 0){
	  cout << "tops have same ID sign!" << endl;
	}
	else{
	  /*
	  cout << "----------------------------------" << endl;
	  cout << "particle " << topID[0] << " pt = " << topPt[0] << endl;
	  cout << "particle " << topID[1] << " pt = " << topPt[1] << endl;
	  */
	  top1.SetPtEtaPhiE(topPt[0],topEta[0],topPhi[0],topEnergy[0]);
	  top2.SetPtEtaPhiE(topPt[1],topEta[1],topPhi[1],topEnergy[1]);
	  
	  ttbarMass = (top1+top2).M();
	}
      }

      // loop over all gen particles in event
      for(size_t i = 0; i < genParticles->size(); i++){
	const reco::GenParticle &p = (*genParticles).at(i);
	int id = p.pdgId();
	
	// find tops 
	if(abs(id) != 6) continue;
	
	allTopsID.push_back(id);
	allTopsPt.push_back(p.pt());
	allTopsEta.push_back(p.eta());
	allTopsPhi.push_back(p.phi());
	allTopsEnergy.push_back(p.energy());
	allTopsStatus.push_back(p.status());
	
      }
    }
    SetValue("topID",topID);
    SetValue("topPt",topPt);
    SetValue("topEta",topEta);
    SetValue("topPhi",topPhi);
    SetValue("topMass",topMass);
    SetValue("topEnergy",topEnergy);
    
    SetValue("ttbarMass",ttbarMass);

    SetValue("allTopsID",allTopsID);
    SetValue("allTopsPt",allTopsPt);
    SetValue("allTopsEta",allTopsEta);
    SetValue("allTopsPhi",allTopsPhi);
    SetValue("allTopsStatus",allTopsStatus);
    SetValue("allTopsEnergy",allTopsEnergy);

    return 0;
}

