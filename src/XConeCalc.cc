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
  bool 	VarNumJets;
  int XConeNumJets_optimal;
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

   // Turn on optimization based on tauDiv 
    if(mPset.exists("VarNumJets")) VarNumJets = mPset.getParameter<bool>("VarNumJets");
    else VarNumJets = true; 

    std::cout << "XConeCalc: VarNumJets = " << VarNumJets << std::endl;
    if(VarNumJets)std::cout << "XConeCalc: IGNORING XConeNumJets ! " << std::endl;

    
   // Define the jet finding plugins for beta = 1.0 , default is 2.0
    if(mPset.exists("XConeBeta")) XConeBeta = mPset.getParameter<double>("XConeBeta");
    else XConeBeta = 2.0; //default

    std::cout << "XConeCalc: XConeBeta = " << XConeBeta << std::endl;

   // Use CHS
    if(mPset.exists("usePFchs")) usePFchs = mPset.getParameter<bool>("usePFchs");
    else usePFchs = true; //default

    std::cout << "XConeCalc: usePFchs = " << usePFchs << std::endl;

   // DEBUG
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
   
   //Things for NJettiness
   double delta;
   if (XConeBeta > 1) delta = 1/(XConeBeta - 1);
   else delta = std::numeric_limits<int>::max(); // use winner take all
   
   double power;
   power = 1.0/XConeBeta;
   
   //NJettiness Stuff - start   
   const int Nmax = 20;
   Njettiness _njettiness(OnePass_GenET_GenKT_Axes(delta, power, XConeR), XConeMeasure(XConeBeta, XConeR));

   //NJettiness nominal
   std::vector<double> tau_njettiness;
   double tau[Nmax];
   for(int N=0; N<Nmax;N++){
	   tau[N] =  _njettiness.getTau(N,FJConstituents);
	   tau_njettiness.push_back( tau[N] );
	   //cout << "njettiness " <<" (N="<< N << ") = " << tau_njettiness.at(N) << endl; 	
   }
   SetValue("tau_njettiness",     tau_njettiness);

   // NJettiness Diff
   std::vector<double> tau_njettiness_diff;
   double tauDiff[Nmax];
   for(int N=0; N<Nmax;N++){
	   if(N==0)tauDiff[N] =  0;
	   if(N!=0)tauDiff[N] =  tau[N]-tau[N-1];
	   tau_njettiness_diff.push_back( tauDiff[N] );
	   //cout << "njettiness_diff " <<" (N="<< N << ") = " << tau_njettiness_diff.at(N) << endl; 	
   }
   SetValue("tau_njettiness_diff",     tau_njettiness_diff);

   // NJettiness Div
   std::vector<double> tau_njettiness_div;
   double tauDiv[Nmax];
   for(int N=0; N<Nmax;N++){
	   if(N==0)tauDiv[N] =  0;
	   if(N!=0)tauDiv[N] =  tau[N]/tau[N-1];
	   tau_njettiness_div.push_back( tauDiv[N] );
	   //cout << "njettiness_div " <<" (N="<< N << ") = " << tau_njettiness_div.at(N) << endl; 	
   }
   SetValue("tau_njettiness_div",     tau_njettiness_div);

   //Finding tauDiv Mean
   Double_t tauDiv_Tot=0;
   Double_t tauDiv_Mean=0;		  
   for(int N=1; N<Nmax;N++){ //start at N=1
	   tauDiv_Tot = tauDiv_Tot+tauDiv[N];
   }
   tauDiv_Mean = tauDiv_Tot /(Nmax-1);
   if(DEBUG)std::cout << "tauDiv_Mean = " << tauDiv_Mean << std::endl;
   
   //Finding N_opt
   int param = 2; //x distance before y below threshold.
   for(int N=Nmax-1; N>-1;N--){ //count from large N
	   if((tauDiv[N]<=tauDiv_Mean)){
		   XConeNumJets_optimal = N + param; //get value just before y is below mean
		   if(DEBUG)std::cout << "XConeNumJets_optimal = " << XConeNumJets_optimal << std::endl;
		   break;
	   }
   }

   //NJettiness Stuff - end   

   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
   if(DEBUG)cout << "Using the XCone Jet Algorithm" << endl;
   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;

   // define the plugins
   
   int N;
   if(VarNumJets) N = XConeNumJets_optimal;
   else N = XConeNumJets;
   
   XConePlugin xcone_pluginA(N, XConeR, XConeBeta);

   // and the jet definitions
   JetDefinition xcone_jetDefA(&xcone_pluginA);

   // and the cluster sequences
   ClusterSequence xcone_seqA(FJConstituents, xcone_jetDefA);

   // and find the jets
   vector<PseudoJet> xcone_jetsA = xcone_seqA.inclusive_jets();

    if(DEBUG)std::cout << "---- " << XConeNumJets <<" XCone Jets ---- R = "<< XConeR << std::endl;
    for (std::vector<fastjet::PseudoJet>::const_iterator ijet = xcone_jetsA.begin(); ijet != xcone_jetsA.end(); ijet++) {
      int index = (int)(ijet-xcone_jetsA.begin());
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
