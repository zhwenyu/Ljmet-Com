/*
 Calculator for substructure variables

 Author: Rizki Syarif 2016. -- preliminary implementation of XCone in LJMet. STILL TESTING PHASE!
 */

#include <iostream>
#include <limits>   // std::numeric_limits
#include <vector>
#include <string>
#include "TLorentzVector.h"

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/BaseEventSelector.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace std;

//Implementing XCone -start
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include <sstream>
#include <ostream>
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh" // In external code, this should be fastjet/contrib/Nsubjettiness.hh
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/XConePlugin.hh"

using namespace fastjet;
using namespace fastjet::contrib;
//Implementing XCone -end


class LjmetFactory;

class XConeCalc : public BaseCalc {

public:
    XConeCalc();
    virtual ~XConeCalc();
    virtual int BeginJob();
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob();

private:
  edm::InputTag packedPFCandColl_it; //added by rizki
  double XConeR;
  int 	XConeNumJets;
  double XConeBeta;
  bool usePFchs;
  bool DEBUG;
	
};

static int reg = LjmetFactory::GetInstance()->Register(new XConeCalc(), "XConeCalc");


XConeCalc::XConeCalc()
{
}

XConeCalc::~XConeCalc()
{
}

int XConeCalc::BeginJob()
{

    if (mPset.exists("packedPFCandColl")){ packedPFCandColl_it = mPset.getParameter<edm::InputTag>("packedPFCandidates");} //added by rizki
    else packedPFCandColl_it = edm::InputTag("packedPFCandidates"); //added by rizki

   // Jet radius to use throughout
    if(mPset.exists("XConeR")) XConeR = mPset.getParameter<double>("XConeR");
    else XConeR = 0.4; 

    std::cout << "XConeCalc: XConeR = " << XConeR << std::endl;


   // Number of Jets to return
    if(mPset.exists("XConeNumJets")) XConeNumJets = mPset.getParameter<int>("XConeNumJets");
    else XConeNumJets = 6; 

    std::cout << "XConeCalc: XConeNumJets = " << XConeNumJets << std::endl;

    
   // Define the jet finding plugins for beta = 1.0 , default is 2.0
    if(mPset.exists("XConeBeta")) XConeBeta = mPset.getParameter<double>("XConeBeta");
    else XConeBeta = 2.0; //default

    std::cout << "XConeCalc: XConeBeta = " << XConeBeta << std::endl;

   // Use CHS
    if(mPset.exists("usePFchs")) usePFchs = mPset.getParameter<bool>("usePFchs");
    else usePFchs = true; //default

    std::cout << "XConeCalc: usePFchs = " << usePFchs << std::endl;

   // DEUG
    if(mPset.exists("DEBUG")) DEBUG = mPset.getParameter<bool>("DEBUG");
    else DEBUG = false; //default

    std::cout << "XConeCalc: DEBUG = " << DEBUG << std::endl;


    return 0;
}

int XConeCalc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector)
{

	//Implementing XCone - start - Rizki
	
	std::vector<double> theXConeJetPt;
    std::vector<double> theXConeJetEta;
    std::vector<double> theXConeJetPhi;
    std::vector<double> theXConeJetEnergy;

	std::vector<fastjet::PseudoJet> FJConstituents;

   edm::Handle<std::vector<pat::PackedCandidate> > PFparticles;
   event.getByLabel(packedPFCandColl_it, PFparticles);
   
   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
   if(DEBUG)cout << "Collecting PFparticles as Jet constituents" << endl;
   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
   if(DEBUG)std::cout << "No. of PFparticles (All) : "<< PFparticles->size() << std::endl;
   int N_PF = 0;
   for (std::vector<pat::PackedCandidate>::const_iterator iPF = PFparticles->begin(); iPF != PFparticles->end(); iPF++) {
      int index = (int)(iPF-PFparticles->begin());
      if(usePFchs){
      	//ATTENTION!! NEED TO CHECK if This is the correct definition for PF CHS!!!
      	if (iPF->fromPV()==0)continue; //CHS - https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#PV_Assignment , https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h#L431
      }
      FJConstituents.push_back( fastjet::PseudoJet( iPF->px(), iPF->py(), iPF->pz(), iPF->energy() ) );
      N_PF++;
    }
   if(DEBUG)std::cout << "No. of PFparticles (as input to XCone): "<< N_PF << std::endl;

   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
   if(DEBUG)cout << "Using the XCone Jet Algorithm" << endl;
   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;

   // define the plugins
   XConePlugin xcone_pluginA_beta1(XConeNumJets, XConeR, XConeBeta);

   // and the jet definitions
   JetDefinition xcone_jetDefA_beta1(&xcone_pluginA_beta1);

   // and the cluster sequences
   ClusterSequence xcone_seqA_beta1(FJConstituents, xcone_jetDefA_beta1);

   // and find the jets
   vector<PseudoJet> xcone_jetsA_beta1 = xcone_seqA_beta1.inclusive_jets();

    if(DEBUG)std::cout << "---- " << XConeNumJets <<" XCone Jets ---- R = "<< XConeR << std::endl;
    for (std::vector<fastjet::PseudoJet>::const_iterator ijet = xcone_jetsA_beta1.begin(); ijet != xcone_jetsA_beta1.end(); ijet++) {
      int index = (int)(ijet-xcone_jetsA_beta1.begin());
      if(DEBUG)std::cout << "no. : " << index << std::endl;
      if(DEBUG)std::cout << "      " <<  ", Jet constituents      : "<< ijet->constituents().size() << std::endl;
      if(DEBUG)std::cout << "      " <<  ", Jet pt                : "<< ijet->pt() << std::endl;
      if(DEBUG)std::cout << "      " <<  ", Jet eta                : "<< ijet->eta() << std::endl;
      if(DEBUG)std::cout << "      " <<  ", Jet phi                : "<< ijet->phi() << std::endl;
      
      theXConeJetPt     . push_back(ijet->pt());
      theXConeJetEta    . push_back(ijet->eta());
      theXConeJetPhi    . push_back(ijet->phi());
      theXConeJetEnergy . push_back(ijet->e());

      }

    SetValue("theXConeJetPt",     theXConeJetPt);
    SetValue("theXConeJetEta",    theXConeJetEta);
    SetValue("theXConeJetPhi",    theXConeJetPhi);
    SetValue("theXConeJetEnergy", theXConeJetEnergy);

	//Implementing XCone - end - Rizki

    return 0;
}

int XConeCalc::EndJob()
{
  return 0;
}
