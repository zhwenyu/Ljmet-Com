/*
  Calculator for topo variables needing Fat Jets

  Author: Joshua Swanson, 2014
*/



#include <iostream>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

class LjmetFactory;

class CATopoCalc : public BaseCalc{
  
public:
  
    CATopoCalc();
    virtual ~CATopoCalc(){}

    virtual int BeginJob(){

        if (mPset.exists("debug"))      debug_ = mPset.getParameter<bool>("debug");
        else                            debug_ = false;

	    if (mPset.exists("CA8PrunedJetColl")) CA8PrunedJetColl_it = mPset.getParameter<edm::InputTag>("CA8PrunedJetColl");
	    else                                  CA8PrunedJetColl_it = edm::InputTag("goodPatJetsCA8PrunedPFPacked");
	
	return 0;
    }
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}

  
private:

    bool debug_;
    edm::InputTag             CA8PrunedJetColl_it;

    int FillBranches( std::vector<edm::Ptr<pat::Muon> > const & vTightMuons,
                           std::vector<edm::Ptr<pat::Electron> > const & vTightElectrons,
                           std::vector<std::pair<TLorentzVector,bool> >  const & vCorrBtagJets,
                           TLorentzVector const & corrMET,
                           std::vector<TLorentzVector> const & vCAWJets,
                           bool isMuon );

};



static int reg = LjmetFactory::GetInstance()->Register(new CATopoCalc(), "CATopoCalc");



CATopoCalc::CATopoCalc(){
  mLegend = "[CATopoCalc]: ";
}


int CATopoCalc::AnalyzeEvent(edm::EventBase const & event,
                                   BaseEventSelector * selector){
    //
    // compute event variables here
    //

    //
    // _____ Get objects from the selector _____________________
    //
    std::vector<edm::Ptr<pat::Jet> >      const & vSelJets = selector->GetSelectedJets();
    std::vector<edm::Ptr<pat::Jet> >      const & vSelBtagJets = selector->GetSelectedBtagJets();
    std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets = selector->GetCorrJetsWithBTags();
    std::vector<edm::Ptr<pat::Jet> >      const & vAllJets = selector->GetAllJets();
    std::vector<edm::Ptr<pat::Muon> >     const & vSelMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons = selector->GetSelectedElectrons();
    edm::Ptr<pat::MET>                    const & pMet = selector->GetMet();
    TLorentzVector                        const & corrMET = selector->GetCorrectedMet();
  

    bool muonchannel = false;
    bool electronchannel = false;
    if (vSelMuons.size()>0) muonchannel = true;
    if (vSelElectrons.size()>0) electronchannel = true;
    if (muonchannel && electronchannel && debug_) 
      std::cout << mLegend <<"WARNING: Two Leptons in the event"<<std::endl;

    if (vSelMuons.size()==0 and vSelElectrons.size()==0) {
      if (debug_) 
	std::cout << mLegend << "No Lepton in event! "<<std::endl;
      return 0;
    }

    if (vCorrBtagJets.size()==0) {
      if (debug_) std::cout << mLegend <<"No Corrected jets in event! "<<std::endl;
        return 0;
    }

    if (corrMET.Pt()==0) {
      if (debug_) std::cout << mLegend <<"No Corrected MET in event! "<<std::endl;
        return 0;
    }

	
    edm::Handle<std::vector<pat::Jet> > CAWJets;
    event.getByLabel(CA8PrunedJetColl_it, CAWJets);
	std::vector<TLorentzVector>  CAWP4;
	TLorentzVector CAJet;
	
    for (std::vector<pat::Jet>::const_iterator ijet = CAWJets->begin(); ijet != CAWJets->end(); ijet++){
    
    	
    	CAJet.SetPxPyPzE( ijet->px(),
    					  ijet->py(),
    					  ijet->pz(),
    					  ijet->energy()	);
    									  
		CAWP4.push_back(CAJet);
	}

    FillBranches(vSelMuons,
                      vSelElectrons,
                      vCorrBtagJets,
                      corrMET,
                      CAWP4,
                      muonchannel, // isMuon                   
                      );
    return 0;
}


int CATopoCalc::FillBranches( std::vector<edm::Ptr<pat::Muon> > const & vSelMuons,
                                         std::vector<edm::Ptr<pat::Electron> > const & vSelElectrons,
                                         std::vector<std::pair<TLorentzVector,bool>> const & vCorrBtagJets, 
                                         TLorentzVector const & corrMET,
                                         std::vector<TLorentzVector> const & vCAWJets,
                                         bool isMuon,                                                    
                                         ){

	TLorentzVector tlv_lepton;

	if ( vSelMuons.size() == 0 && vSelElectrons.size() == 0) break;
	if ( vSelMuons.size() > 0 && vSelElectrons.size() > 0) {
		tlv_lepton.SetPxPyPzE( vSelMuons[0]->px(),
							   vSelMuons[0]->py(),
							   vSelMuons[0]->pz(),
							   vSelMuons[0]->energy() );
	}
	if ( vSelMuons.size() > 0 ) {
		tlv_lepton.SetPxPyPzE( vSelMuons[0]->px(),
							   vSelMuons[0]->py(),
							   vSelMuons[0]->pz(),
							   vSelMuons[0]->energy() ); 
	}
	if ( vSelElectrons.size() > 0 ) {
		tlv_lepton.SetPxPyPzE( vSelElectrons[0]->px(),
							   vSelElectrons[0]->py(),
							   vSelElectrons[0]->pz(),
							   vSelElectrons[0]->energy() );
	}
	
	if ( corrMET.Pt() > 0 ){ }
	else break;

	TLorentzVector tlv_met( corrMET.Px(),
							corrMET.Py(),
							corrMET.Pz(),
							corrMET.Energy() );

	std::vector<TLorentzVector> jets;
	std::vector<TLorentzVector> bjets;
	 
	for (vector<std::pair<TLorentzVector,bool>>::const_iterator jet = vCorrBtagJets.begin(); jet != vCorrBtagJets.end(); ++jet){		
	
		double CAtoAKJetDR = deltaR(jet.first(),vCAWJets[0]);		
		if( CAtoAKJetDR > 0.8 ){		
			if(jet.second)	bjets.push_back(jet);
			else	jets.push_back(jet);	
		}
	}
	double tPrimeMass = -10.;
	if( bjets.size() > 0 && vCAWJets.size() > 0 ){
		tPrimeMass = ( tlv_met + tlv_lepton + vCAWJets[0] + bjets[0] ).M;
	}
	SetValue("tPrimeMass", tPrimeMass);

	return result;


}
