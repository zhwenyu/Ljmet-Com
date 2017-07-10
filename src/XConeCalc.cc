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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

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

#include "DataFormats/Math/interface/deltaR.h"

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
  edm::InputTag packedGenParticleColl_it; //added by rizki
  double	XConeR;
  int		XConeNumJets;
  bool		VarNumJets;
  double	XConeBeta;
  bool		usePFchs;
  bool		doPUPPI;
  bool		doGenXCone;
  bool		saveJetConst;
  bool		isMc;
  bool		DEBUG;
  
  bool		saveLooseLeps;
	
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

    if (mPset.exists("packedGenParticleColl")){ packedGenParticleColl_it = mPset.getParameter<edm::InputTag>("packedGenParticles");} //added by rizki
    else packedGenParticleColl_it = edm::InputTag("packedGenParticles"); //added by rizki

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

   // do XCone PUPPI
    if(mPset.exists("doPUPPI")) doPUPPI = mPset.getParameter<bool>("doPUPPI");
    else doPUPPI = true; //default
    std::cout << "XConeCalc: doXConePUPPI = " << doPUPPI << std::endl;

   // do XCone Gen
    if(mPset.exists("doGenXCone")) doGenXCone = mPset.getParameter<bool>("doGenXCone");
    else doGenXCone = true; //default
    std::cout << "XConeCalc: doXConeGen = " << doGenXCone << std::endl;

   // isMc?
    if(mPset.exists("isMc")) isMc = mPset.getParameter<bool>("isMc");
    else isMc = false;

   // save Jet constituents?
    if(mPset.exists("saveJetConst")) saveJetConst = mPset.getParameter<bool>("saveJetConst");
    else saveJetConst = false; //default
    std::cout << "XConeCalc: saveJetConst = " << saveJetConst << std::endl;

   // DEBUG
    if(mPset.exists("DEBUG")) DEBUG = mPset.getParameter<bool>("DEBUG");
    else DEBUG = false; //default
    std::cout << "XConeCalc: DEBUG = " << DEBUG << std::endl;

    if (mPset.exists("saveLooseLeps")) saveLooseLeps = mPset.getParameter<bool>("saveLooseLeps");
    else                               saveLooseLeps = false;

    return 0;
}

int XConeCalc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector)
{
    std::vector<edm::Ptr<pat::Muon> >           const & vSelectedMuons = selector->GetSelectedMuons();
    std::vector<edm::Ptr<pat::Electron> >       const & vSelectedElectrons = selector->GetSelectedElectrons();
    std::vector<edm::Ptr<pat::Muon> >           const & vLooseMuons = selector->GetLooseMuons();
    std::vector<edm::Ptr<pat::Electron> >       const & vLooseElectrons = selector->GetLooseElectrons();

    std::vector<edm::Ptr<pat::Muon> > vSelMuons;
    std::vector<edm::Ptr<pat::Electron> > vSelElectrons;
    if(saveLooseLeps){
      vSelMuons = vLooseMuons;
      vSelElectrons = vLooseElectrons;
    }else{
      vSelMuons = vSelectedMuons;
      vSelElectrons = vSelectedElectrons;
    }


	//Implementing XCone - start - Rizki
	
	std::vector<fastjet::PseudoJet> FJConstituents; 
  	int		XConeNumJets_optimal=0;
	std::vector<double> theXConeJetPt;
    std::vector<double> theXConeJetEta;
    std::vector<double> theXConeJetPhi;
    std::vector<double> theXConeJetEnergy;
    std::vector<double> theXConeJetArea;
    std::vector<double> theXConeJetConstStartIndex;
    std::vector<double> theXConeJetConstEndIndex;
    //collect constituent info - start
	std::vector<double> theXConeJetConstPt;
    std::vector<double> theXConeJetConstEta;
	std::vector<double> theXConeJetConstPhi;
	std::vector<double> theXConeJetConstEnergy;
    //collect constituent info - end

	//PUPPI:
	std::vector<fastjet::PseudoJet> FJConstituentsPUPPI;
  	int		XConeNumJets_optimal_puppi=0;
	std::vector<double> theXConePUPPIJetPt;
    std::vector<double> theXConePUPPIJetEta;
    std::vector<double> theXConePUPPIJetPhi;
    std::vector<double> theXConePUPPIJetEnergy;
    std::vector<double> theXConePUPPIJetArea;
    std::vector<double> theXConePUPPIJetConstStartIndex;
    std::vector<double> theXConePUPPIJetConstEndIndex;
    //collect constituent info - start
	std::vector<double> theXConePUPPIJetConstPt;
    std::vector<double> theXConePUPPIJetConstEta;
	std::vector<double> theXConePUPPIJetConstPhi;
	std::vector<double> theXConePUPPIJetConstEnergy;
    //collect constituent info - end

	//GenXConeJet:
	std::vector<fastjet::PseudoJet> FJConstituentsGen;
	std::vector<double> theXConeGenJetPt;
    std::vector<double> theXConeGenJetEta;
    std::vector<double> theXConeGenJetPhi;
    std::vector<double> theXConeGenJetEnergy;
    std::vector<double> theXConeGenJetArea;
    std::vector<double> theXConeGenJetConstStartIndex;
    std::vector<double> theXConeGenJetConstEndIndex;
    //collect constituent info - start
	std::vector<double> theXConeGenJetConstPt;
    std::vector<double> theXConeGenJetConstEta;
	std::vector<double> theXConeGenJetConstPhi;
	std::vector<double> theXConeGenJetConstEnergy;
    //collect constituent info - end

	//PF particles
   edm::Handle<std::vector<pat::PackedCandidate> > PFparticles;
   event.getByLabel(packedPFCandColl_it, PFparticles);
   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
   if(DEBUG)cout << "Collecting PFparticles as Jet constituents" << endl;
   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
   if(DEBUG)std::cout << "No. of PFparticles (All) : "<< PFparticles->size() << std::endl;
   int N_PF = 0;
   int N_PFch = 0;
   int N_PFpuppi = 0;
   for (std::vector<pat::PackedCandidate>::const_iterator iPF = PFparticles->begin(); iPF != PFparticles->end(); iPF++) {
      int index = (int)(iPF-PFparticles->begin());
      
//       TLorentzVector iPF_lv;
//       iPF_lv.SetEtaPhiE(iPF->px(), iPF->py(), iPF->pz(), iPF->energy())
      
      //attempt to exclude selectedLeptons PF candidates - start
      bool isPFlep=false;
      double minLepPF_dR = 10000.;
      double lepPF_dR = 0.;
      int i_mu = 0;
      for (std::vector<edm::Ptr<pat::Muon> >::const_iterator imu = vSelMuons.begin(); imu != vSelMuons.end(); imu++) {
//       	double lepM = 0.105658367;
//       	TLorentzVector lep_lv;
//       	lep_lv.SetPtEtaPhiM((*imu)->pt(),(*imu)->eta(),(*imu)->phi(),lepM);
		lepPF_dR = deltaR((*imu)->eta(),(*imu)->phi(),iPF->eta(),iPF->phi());
      	if(lepPF_dR<minLepPF_dR) minLepPF_dR = lepPF_dR;      		
      	if(lepPF_dR < 0.01){
      		if(DEBUG)std::cout << "iPF : "<< index << " ( pT = "<< iPF->pt() << ", eta = "<< iPF->eta() <<", phi = "<< iPF->phi() <<", E = "<< iPF->energy() <<" )" ; 
      		if(DEBUG)std::cout << "		i_mu : "<< i_mu <<"(pT = "<< (*imu)->pt() <<", eta = "<<  (*imu)->eta()<<", phi = "<< (*imu)->phi() <<", E = "<< (*imu)->energy() <<"), minMuPF_dR = "<< minLepPF_dR << ", PF inv mass = " << iPF->p4().M() << std::endl;
      		isPFlep = true;
      	} 
      	i_mu++;
      }
      if(isPFlep){
      	if(DEBUG) cout << "			---> Found lepton match. not clustering into XCone jet!" << endl;
      	continue;
      }		
      int i_el = 0;
      for (std::vector<edm::Ptr<pat::Electron> >::const_iterator iel = vSelElectrons.begin(); iel != vSelElectrons.end(); iel++){
// 		double lepM = 0.00051099891;
//       	TLorentzVector lep_lv;
//       	lep_lv.SetPtEtaPhiM((*iel)->pt(),(*iel)->eta(),(*iel)->phi(),lepM);
		lepPF_dR = deltaR((*iel)->eta(),(*iel)->phi(),iPF->eta(),iPF->phi());
      	if(lepPF_dR<minLepPF_dR) minLepPF_dR = lepPF_dR;   
      	if(lepPF_dR < 0.01){
      		if(DEBUG)std::cout << "iPF : "<< index << " ( pT = "<< iPF->pt() << ", eta = "<< iPF->eta() <<", phi = "<< iPF->phi() <<", E = "<< iPF->energy() <<" )" ; 
      		if(DEBUG)std::cout << "		i_el : "<< i_el <<"(pT = "<< (*iel)->pt() <<", eta = "<<  (*iel)->eta()<<", phi = "<< (*iel)->phi() <<", E = "<< (*iel)->energy() <<"), minMuPF_dR = "<< minLepPF_dR << ", PF inv mass = " << iPF->p4().M() << std::endl;
      		isPFlep = true;
      	} 
      	i_el++;
      }
      if(isPFlep){
      	if(DEBUG) cout << "			---> Found lepton match. not clustering into XCone jet!" << endl;
      	continue;
      }		
      //attempt to exclude selectedLeptons PF candidates - end
      N_PF++;
      
      if(doPUPPI){
		  float wPup = iPF->puppiWeight();
		  //if(DEBUG) cout << "PUPPI weight for PF no." << index << ": "<< wPup << endl;
		  if(wPup!=0){
			  FJConstituentsPUPPI.push_back( fastjet::PseudoJet( iPF->px()*wPup, iPF->py()*wPup, iPF->pz()*wPup, iPF->energy()*wPup ) );
			  N_PFpuppi++;
		  }
		}

      if(usePFchs){
      	//ATTENTION!! NEED TO CHECK if This is the correct definition for PF CHS!!!
      	if (iPF->fromPV()==0){
      		N_PFch++;
      		continue; //CHS - https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#PV_Assignment , https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/PackedCandidate.h#L431
      	}
      }      
      FJConstituents.push_back( fastjet::PseudoJet( iPF->px(), iPF->py(), iPF->pz(), iPF->energy() ) );

    }
   if(DEBUG)std::cout << "No. of PFparticles (after lepton cleaning):	"<< N_PF << std::endl;
   if(DEBUG)std::cout << "No. of PFPUPPIparticles (after lepton cleaning):	"<< N_PFpuppi << std::endl;
   if(DEBUG)std::cout << "No. of PFparticles (after CHS):	"<< N_PF-N_PFch << std::endl;
   
   //GEN PArticles
   edm::Handle<std::vector<pat::PackedGenParticle> > Genparticles;
   event.getByLabel(packedGenParticleColl_it, Genparticles);   
   if(isMc){
	   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
	   if(DEBUG)cout << "Collecting Genparticles as Jet constituents (For MC)" << endl;
	   if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
	   if(DEBUG)std::cout << "No. of Genparticles (All) : "<< Genparticles->size() << std::endl;
	   int N_Gen = 0;
	   for (std::vector<pat::PackedGenParticle>::const_iterator iGen = Genparticles->begin(); iGen != Genparticles->end(); iGen++) {
		  int index = (int)(iGen-Genparticles->begin());
			
		  //attempt to exclude selectedLeptons Gen candidates - start
		  bool isGenlep=false;
		  double minLepGen_dR = 10000.;
		  double lepGen_dR = 0.;
		  int i_mu = 0;
		  for (std::vector<edm::Ptr<pat::Muon> >::const_iterator imu = vSelMuons.begin(); imu != vSelMuons.end(); imu++) {
			if(abs(iGen->pdgId())!=11 && abs(iGen->pdgId())!=13) continue; //only check if Gen is mu/el
			lepGen_dR = deltaR((*imu)->eta(),(*imu)->phi(),iGen->eta(),iGen->phi());
			if(lepGen_dR<minLepGen_dR) minLepGen_dR = lepGen_dR;      		
			if(lepGen_dR < 0.01){
				if(DEBUG)std::cout << "iGen : "<< index << " ( pT = "<< iGen->pt() << ", eta = "<< iGen->eta() <<", phi = "<< iGen->phi() <<", E = "<< iGen->energy() <<" pdgId =" << iGen->pdgId() <<" )" ; 
				if(DEBUG)std::cout << "		i_mu : "<< i_mu <<"(pT = "<< (*imu)->pt() <<", eta = "<<  (*imu)->eta()<<", phi = "<< (*imu)->phi() <<", E = "<< (*imu)->energy() <<"), minMuGen_dR = "<< minLepGen_dR << ", Gen inv mass = " << iGen->p4().M() << std::endl;
				isGenlep = true;
			} 
			i_mu++;
		  }
		  if(isGenlep){
			if(DEBUG) cout << "			---> Found lepton match. not clustering into XCone jet!" << endl;
			continue;
		  }		
		  int i_el = 0;
		  for (std::vector<edm::Ptr<pat::Electron> >::const_iterator iel = vSelElectrons.begin(); iel != vSelElectrons.end(); iel++){
			if(abs(iGen->pdgId())!=11 && abs(iGen->pdgId())!=13) continue; //only check if Gen is mu/el
			lepGen_dR = deltaR((*iel)->eta(),(*iel)->phi(),iGen->eta(),iGen->phi());
			if(lepGen_dR<minLepGen_dR) minLepGen_dR = lepGen_dR;   
			if(lepGen_dR < 0.01){
				if(DEBUG)std::cout << "iGen : "<< index << " ( pT = "<< iGen->pt() << ", eta = "<< iGen->eta() <<", phi = "<< iGen->phi() <<", E = "<< iGen->energy() <<" pdgId =" << iGen->pdgId() <<" )" ; 
				if(DEBUG)std::cout << "		i_el : "<< i_el <<"(pT = "<< (*iel)->pt() <<", eta = "<<  (*iel)->eta()<<", phi = "<< (*iel)->phi() <<", E = "<< (*iel)->energy() <<"), minMuGen_dR = "<< minLepGen_dR << ", Gen inv mass = " << iGen->p4().M() << std::endl;
				isGenlep = true;
			} 
			i_el++;
		  }
		  if(isGenlep){
			if(DEBUG) cout << "			---> Found lepton match. not clustering into XCone jet!" << endl;
			continue;
		  }		
		  //attempt to exclude selectedLeptons Gen candidates - end
		  N_Gen++;
	  
		  FJConstituentsGen.push_back( fastjet::PseudoJet( iGen->px(), iGen->py(), iGen->pz(), iGen->energy() ) );

		}
	   if(DEBUG)std::cout << "No. of Genparticles (after lepton matching):	"<< N_Gen << std::endl;
   }

   
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
   std::vector<double> tau_njettiness_puppi;
   double tau[Nmax];
   double tau_puppi[Nmax];
   for(int N=0; N<Nmax;N++){
	   tau[N] =  _njettiness.getTau(N,FJConstituents);
	   tau_njettiness.push_back( tau[N] );
	   //if(DEBUG)cout << "njettiness " <<" (N="<< N << ") = " << tau_njettiness.at(N) << endl; 	

	   //PUPPI
	   if(!doPUPPI) continue;
	   tau_puppi[N] =  _njettiness.getTau(N,FJConstituentsPUPPI);
	   tau_njettiness_puppi.push_back( tau_puppi[N] );
	   //if(DEBUG)cout << "njettiness_puppi " <<" (N="<< N << ") = " << tau_njettiness_puppi.at(N) << endl; 	
   }
   SetValue("tau_njettiness",     tau_njettiness);
   if(doPUPPI)SetValue("tau_njettiness_puppi",     tau_njettiness_puppi);

   // NJettiness Diff
   std::vector<double> tau_njettiness_diff;
   std::vector<double> tau_njettiness_diff_puppi;
   double tauDiff[Nmax];
   double tauDiff_puppi[Nmax];
   for(int N=0; N<Nmax;N++){
	   if(N==0)tauDiff[N] =  0;
	   if(N!=0)tauDiff[N] =  tau[N]-tau[N-1];
	   tau_njettiness_diff.push_back( tauDiff[N] );
	   //if(DEBUG)cout << "njettiness_diff " <<" (N="<< N << ") = " << tau_njettiness_diff.at(N) << endl; 	

	   //PUPPI
	   if(!doPUPPI) continue;
	   if(N==0)tauDiff_puppi[N] =  0;
	   if(N!=0)tauDiff_puppi[N] =  tau_puppi[N]-tau_puppi[N-1];
	   tau_njettiness_diff_puppi.push_back( tauDiff_puppi[N] );
	   //if(DEBUG)cout << "njettiness_diff_puppi " <<" (N="<< N << ") = " << tau_njettiness_diff_puppi.at(N) << endl; 	
   }
   SetValue("tau_njettiness_diff",     tau_njettiness_diff);
   if(doPUPPI)SetValue("tau_njettiness_diff_puppi",     tau_njettiness_diff_puppi);

   // NJettiness Div
   std::vector<double> tau_njettiness_div;
   std::vector<double> tau_njettiness_div_puppi;
   double tauDiv[Nmax];
   double tauDiv_puppi[Nmax];
   for(int N=0; N<Nmax;N++){
	   if(N==0)tauDiv[N] =  0;
	   if(N!=0)tauDiv[N] =  tau[N]/tau[N-1];
	   tau_njettiness_div.push_back( tauDiv[N] );
	   //if(DEBUG)cout << "njettiness_div " <<" (N="<< N << ") = " << tau_njettiness_div.at(N) << endl; 	

	   //PUPPI
	   if(!doPUPPI) continue;
	   if(N==0)tauDiv_puppi[N] =  0;
	   if(N!=0)tauDiv_puppi[N] =  tau_puppi[N]/tau_puppi[N-1];
	   tau_njettiness_div_puppi.push_back( tauDiv_puppi[N] );
	   //if(DEBUG)cout << "njettiness_div_puppi " <<" (N="<< N << ") = " << tau_njettiness_div_puppi.at(N) << endl; 	
   }
   SetValue("tau_njettiness_div",     tau_njettiness_div);
   if(doPUPPI)SetValue("tau_njettiness_div_puppi",     tau_njettiness_div_puppi);

   // NJettiness DiffDivComb
   std::vector<double> tau_njettiness_diffdivComb;
   std::vector<double> tau_njettiness_diffdivComb_puppi;
   double tauDiffDivComb[Nmax];
   double tauDiffDivComb_puppi[Nmax];
   for(int N=0; N<Nmax;N++){
	   if(N==0)tauDiffDivComb[N] =  0;
	   if(N!=0)tauDiffDivComb[N] =  (tau[N]-tau[N-1]) /tau[N];
	   tau_njettiness_diffdivComb.push_back( tauDiffDivComb[N] );
	   //if(DEBUG)cout << "njettiness_diffdivComb " <<" (N="<< N << ") = " << tau_njettiness_diffdivComb.at(N) << endl;
	   
	   //PUPPI
	   if(!doPUPPI) continue;
	   if(N==0)tauDiffDivComb_puppi[N] =  0;
	   if(N!=0)tauDiffDivComb_puppi[N] =  (tau_puppi[N]-tau_puppi[N-1]) /tau_puppi[N];
	   tau_njettiness_diffdivComb_puppi.push_back( tauDiffDivComb_puppi[N] );
	   //if(DEBUG)cout << "njettiness_diffdivComb_puppi " <<" (N="<< N << ") = " << tau_njettiness_diffdivComb_puppi.at(N) << endl; 	
   }
   SetValue("tau_njettiness_diffdivComb",     tau_njettiness_diffdivComb);
   if(doPUPPI)SetValue("tau_njettiness_diffdivComb_puppi",     tau_njettiness_diffdivComb_puppi);
   
   /*
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
   */

   //Finding N_opt
   if(DEBUG) cout << "Finding N_opt ... " << endl;
   double f_opt[Nmax]; for(int N=0;N<Nmax;N++)f_opt[N]=tauDiff[N]; double f_max = 0. ; if(DEBUG) cout << " using tauDiff " << endl;
   double reduceThresh = 5; //need to be configurable!
   double Thresh = f_max-30.; //need to be configurable!
//    double f_opt[Nmax]; for(int N=0;N<Nmax;N++)f_opt[N]=tauDiv[N]; double f_max = 1. ; if(DEBUG) cout << " using tauDiv " << endl;
//    double f_opt[Nmax]; for(int N=0;N<Nmax;N++)f_opt[N]= tauDiffDivComb; double f_max = 0. ; if(DEBUG) cout << " using tauDiffDivComb " << endl;
//    double reduceThresh = 0.025; //need to be configurable!
//    double Thresh = 0.025; //need to be configurable!
   if(DEBUG) cout <<"Threshold = "<< Thresh << endl;
   int param = 0; //distance after N pass.
   XConeNumJets_optimal = 0.;
   while(XConeNumJets_optimal==0.){
	   for(int N=1; N<Nmax; N++){
		   if(DEBUG) cout << "f_opt["<<N<<"] = " << f_opt[N] << endl;
		   if((f_opt[N]>Thresh)){
			   XConeNumJets_optimal = N + param; 
			   break;
		   }
	   }
	   Thresh = Thresh - reduceThresh; //lower theshold if N_opt is not found
	   if(DEBUG) cout << "XConeNumJets_optimal = " << XConeNumJets_optimal << endl;
   }
   SetValue("XConeNumJets_optimal",     XConeNumJets_optimal);


   //Finding N_opt with PUPPI
   if(doPUPPI){
	   if(DEBUG) cout << "Finding N_opt (PUPPI) ... " << endl;
	   double f_opt_puppi[Nmax]; for(int N=0;N<Nmax;N++)f_opt_puppi[N]=tauDiff_puppi[N]; double f_max_puppi = 0. ; if(DEBUG) cout << " using tauDiff (for PUPPI)" << endl;
// 	   double f_opt_puppi[Nmax]; for(int N=0;N<Nmax;N++)f_opt_puppi[N]=tauDiv_puppi[N]; double f_max_puppi = 1. ; if(DEBUG) cout << " using tauDiv (for PUPPI)" << endl;
// 	   double f_opt_puppi[Nmax]; for(int N=0;N<Nmax;N++)f_opt_puppi[N]= tauDiffDivComb_puppi; double f_max_puppi = 0. ; if(DEBUG) cout << " using tauDiffDivComb (for PUPPI) " << endl;
	   double reduceThresh_puppi = 5.; //need to be configurable!
	   double Thresh_puppi = -25.; //need to be configurable!
	   if(DEBUG) cout <<"Threshold_puppi = "<< Thresh_puppi << endl;
	   int param_puppi = 0; //distance after N pass.
	   XConeNumJets_optimal_puppi = 0.;
	   while(XConeNumJets_optimal_puppi==0.){
		   for(int N=1; N<Nmax; N++){
			   if(DEBUG) cout << "f_opt_puppi["<<N<<"] = " << f_opt_puppi[N] << endl;
			   if((f_opt_puppi[N]>Thresh_puppi)){
				   XConeNumJets_optimal_puppi = N + param_puppi; 
				   break;
			   }
		   }
		   Thresh_puppi = Thresh_puppi - reduceThresh_puppi; //lower theshold if N_opt is not found
		   if(DEBUG) cout << "XConeNumJets_optimal_puppi = " << XConeNumJets_optimal_puppi << endl;
	   }
	   SetValue("XConeNumJets_optimal_puppi",     XConeNumJets_optimal_puppi);
   }

   //NJettiness Stuff - end   

	//Attempting to to define jet area - start
	GhostedAreaSpec ghost_spec(3); //set max rapiditiy for creating ghosts.
	AreaDefinition areaDef(active_area_explicit_ghosts,ghost_spec);   
	//Attempting to to define jet area - end

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
//    ClusterSequence xcone_seqA(FJConstituents, xcone_jetDefA);
   ClusterSequenceArea xcone_seqA(FJConstituents, xcone_jetDefA, areaDef);

   // and find the jets
   vector<PseudoJet> xcone_jetsA = xcone_seqA.inclusive_jets();

    if(DEBUG)std::cout << "---- " << N <<" XCone Jets ---- R = "<< XConeR << std::endl;
    int ConstituentStartIndex =0;
    for (std::vector<fastjet::PseudoJet>::const_iterator ijet = xcone_jetsA.begin(); ijet != xcone_jetsA.end(); ijet++) {
      int index = (int)(ijet-xcone_jetsA.begin());
      if(DEBUG)std::cout << "no. : " << index << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet constituents	: "<< ijet->constituents().size() << "	("<<ConstituentStartIndex <<" - "<< ConstituentStartIndex+ijet->constituents().size()-1 << ") "<<std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet pt			: "<< ijet->pt() << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet eta			: "<< ijet->eta() << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet phi			: "<< ijet->phi() << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet area		: "<< ijet->area() << std::endl;
      
      theXConeJetPt     . push_back(ijet->pt());
      theXConeJetEta    . push_back(ijet->eta());
      theXConeJetPhi    . push_back(ijet->phi());
      theXConeJetEnergy . push_back(ijet->e());
      theXConeJetArea   . push_back(ijet->area());
      
      //collect constituent info - start
      if(saveJetConst){
		  theXConeJetConstStartIndex     . push_back(ConstituentStartIndex);
		  theXConeJetConstEndIndex     . push_back(ConstituentStartIndex+ijet->constituents().size()-1);
		  ConstituentStartIndex+=ijet->constituents().size();

		  for (unsigned int iconst = 0; iconst < ijet->constituents().size(); iconst++) {
				theXConeJetConstPt     . push_back(ijet->constituents().at(iconst).pt());
				theXConeJetConstEta    . push_back(ijet->constituents().at(iconst).eta());
				theXConeJetConstPhi    . push_back(ijet->constituents().at(iconst).phi());
				theXConeJetConstEnergy . push_back(ijet->constituents().at(iconst).e());
		  }
		  /*if(DEBUG){
			for(unsigned int i=theXConeJetConstStartIndex.at(index); i<= theXConeJetConstEndIndex.at(index);i++){
				std::cout << "		const no: "<< i ;
				std::cout << "	pt: "<<theXConeJetConstPt.at(i);
				std::cout << "	eta: "<<theXConeJetConstEta.at(i);
				std::cout << "	phi: "<<theXConeJetConstPhi.at(i);
				std::cout << "	energy: "<<theXConeJetConstEnergy.at(i); 
				if(theXConeJetConstPt.at(i)<1e-50)std::cout << "	-----> GHOST!!"; 
				std::cout << std::endl;
			}
		  }*/
	  }	
      //collect constituent info - end

    }

    SetValue("theXConeJetPt",     theXConeJetPt);
    SetValue("theXConeJetEta",    theXConeJetEta);
    SetValue("theXConeJetPhi",    theXConeJetPhi);
    SetValue("theXConeJetEnergy", theXConeJetEnergy);
    SetValue("theXConeJetArea", theXConeJetArea);

    SetValue("theXConeJetConstStartIndex",     theXConeJetConstStartIndex);
    SetValue("theXConeJetConstEndIndex",    theXConeJetConstEndIndex);
    SetValue("theXConeJetConstPt",     theXConeJetConstPt);
    SetValue("theXConeJetConstEta",    theXConeJetConstEta);
    SetValue("theXConeJetConstPhi",    theXConeJetConstPhi);
    SetValue("theXConeJetConstEnergy", theXConeJetConstEnergy);

	if(doPUPPI){

		if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
		if(DEBUG)cout << "Using the XCone Jet Algorithm (PUPPI) " << endl;
		if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;

		// define the plugins

		int N_puppi;
		if(VarNumJets) N_puppi = XConeNumJets_optimal_puppi;
		else N_puppi = XConeNumJets;

		XConePlugin xcone_pluginA_puppi(N_puppi, XConeR, XConeBeta);

		// and the jet definitions
		JetDefinition xcone_jetDefA_puppi(&xcone_pluginA_puppi);

		// and the cluster sequences
		// ClusterSequence xcone_seqA_puppi(FJConstituents_puppi, xcone_jetDefA_puppi);
		ClusterSequenceArea xcone_seqA_puppi(FJConstituentsPUPPI, xcone_jetDefA_puppi, areaDef);

		// and find the jets (PUPPI)
		vector<PseudoJet> xcone_jetsA_puppi = xcone_seqA_puppi.inclusive_jets();

		if(DEBUG)std::cout << "---- " << N_puppi <<" XConePUPPI Jets ---- R = "<< XConeR << std::endl;
		int PUPPIConstituentStartIndex = 0;
		for (std::vector<fastjet::PseudoJet>::const_iterator ijet = xcone_jetsA_puppi.begin(); ijet != xcone_jetsA_puppi.end(); ijet++) {
		  int index = (int)(ijet-xcone_jetsA_puppi.begin());
		  if(DEBUG)std::cout << "no. : " << index << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet constituents	: "<< ijet->constituents().size() << "	("<<PUPPIConstituentStartIndex <<" - "<< PUPPIConstituentStartIndex+ijet->constituents().size()-1 << ") "<<std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet pt			: "<< ijet->pt() << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet eta			: "<< ijet->eta() << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet phi			: "<< ijet->phi() << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet area		: "<< ijet->area() << std::endl;
  
		  theXConePUPPIJetPt     . push_back(ijet->pt());
		  theXConePUPPIJetEta    . push_back(ijet->eta());
		  theXConePUPPIJetPhi    . push_back(ijet->phi());
		  theXConePUPPIJetEnergy . push_back(ijet->e());
		  theXConePUPPIJetArea   . push_back(ijet->area());

		  //collect constituent info - start
		  if(saveJetConst){
			  theXConePUPPIJetConstStartIndex     . push_back(PUPPIConstituentStartIndex);
			  theXConePUPPIJetConstEndIndex     . push_back(PUPPIConstituentStartIndex+ijet->constituents().size()-1);
			  PUPPIConstituentStartIndex+=ijet->constituents().size();
			  for (unsigned int iconst = 0; iconst < ijet->constituents().size(); iconst++) {
					theXConePUPPIJetConstPt     . push_back(ijet->constituents().at(iconst).pt());
					theXConePUPPIJetConstEta    . push_back(ijet->constituents().at(iconst).eta());
					theXConePUPPIJetConstPhi    . push_back(ijet->constituents().at(iconst).phi());
					theXConePUPPIJetConstEnergy . push_back(ijet->constituents().at(iconst).e());
			  }
			  /*if(DEBUG){
				for(unsigned int i=theXConePUPPIJetConstStartIndex.at(index); i<= theXConePUPPIJetConstEndIndex.at(index);i++){
					std::cout << "		const no: "<< i ;
					std::cout << "	pt: "<<theXConePUPPIJetConstPt.at(i);
					std::cout << "	eta: "<<theXConePUPPIJetConstEta.at(i);
					std::cout << "	phi: "<<theXConePUPPIJetConstPhi.at(i);
					std::cout << "	energy: "<<theXConePUPPIJetConstEnergy.at(i); 
					if(theXConePUPPIJetConstPt.at(i)<1e-50)std::cout << "	-----> GHOST!!"; 
					std::cout << std::endl;
				}
			  }*/
		  }
		  //collect constituent info - end

		}

		SetValue("theXConePUPPIJetPt",     theXConePUPPIJetPt);
		SetValue("theXConePUPPIJetEta",    theXConePUPPIJetEta);
		SetValue("theXConePUPPIJetPhi",    theXConePUPPIJetPhi);
		SetValue("theXConePUPPIJetEnergy", theXConePUPPIJetEnergy);
		SetValue("theXConePUPPIJetArea", theXConePUPPIJetArea);

		SetValue("theXConePUPPIJetConstStartIndex",     theXConePUPPIJetConstStartIndex);
		SetValue("theXConePUPPIJetConstEndIndex",     theXConePUPPIJetConstEndIndex);
		SetValue("theXConePUPPIJetConstPt",     theXConePUPPIJetConstPt);
		SetValue("theXConePUPPIJetConstEta",    theXConePUPPIJetConstEta);
		SetValue("theXConePUPPIJetConstPhi",    theXConePUPPIJetConstPhi);
		SetValue("theXConePUPPIJetConstEnergy", theXConePUPPIJetConstEnergy);

	}

	if(isMc & doGenXCone){

		if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;
		if(DEBUG)cout << "Using the XCone Jet Algorithm (Gen) " << endl;
		if(DEBUG)cout << "-------------------------------------------------------------------------------------" << endl;

		// define the plugins

		int N_Gen;
		if(VarNumJets){
			N_Gen = XConeNumJets_optimal; //Sync with PFchs N, is this ok?? Need to come back to this and check!
			if(DEBUG)cout << "using optimal N from PFchs!" << endl;
		}
		else N_Gen = XConeNumJets;

		XConePlugin xcone_pluginA_Gen(N_Gen, XConeR, XConeBeta);

		// and the jet definitions
		JetDefinition xcone_jetDefA_Gen(&xcone_pluginA_Gen);

		// and the cluster sequences
		// ClusterSequence xcone_seqA_Gen(FJConstituents_Gen, xcone_jetDefA_Gen);
		ClusterSequenceArea xcone_seqA_Gen(FJConstituentsGen, xcone_jetDefA_Gen, areaDef);

		// and find the jets (Gen)
		vector<PseudoJet> xcone_jetsA_Gen = xcone_seqA_Gen.inclusive_jets();

		if(DEBUG)std::cout << "---- " << N_Gen <<" XConeGen Jets ---- R = "<< XConeR << std::endl;
    	int GenConstituentStartIndex =0;
		for (std::vector<fastjet::PseudoJet>::const_iterator ijet = xcone_jetsA_Gen.begin(); ijet != xcone_jetsA_Gen.end(); ijet++) {
		  int index = (int)(ijet-xcone_jetsA_Gen.begin());
		  if(DEBUG)std::cout << "no. : " << index << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet constituents	: "<< ijet->constituents().size() << "	("<<GenConstituentStartIndex <<" - "<< GenConstituentStartIndex+ijet->constituents().size()-1 << ") "<<std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet pt			: "<< ijet->pt() << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet eta			: "<< ijet->eta() << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet phi			: "<< ijet->phi() << std::endl;
		  if(DEBUG)std::cout << "      " <<  "Jet area		: "<< ijet->area() << std::endl;
  
		  theXConeGenJetPt     . push_back(ijet->pt());
		  theXConeGenJetEta    . push_back(ijet->eta());
		  theXConeGenJetPhi    . push_back(ijet->phi());
		  theXConeGenJetEnergy . push_back(ijet->e());
		  theXConeGenJetArea   . push_back(ijet->area());

		  //collect constituent info - start
		  if(saveJetConst){		  
			  theXConeGenJetConstStartIndex     . push_back(GenConstituentStartIndex);
			  theXConeGenJetConstEndIndex     . push_back(GenConstituentStartIndex+ijet->constituents().size()-1);
			  GenConstituentStartIndex+=ijet->constituents().size();
			  for (unsigned int iconst = 0; iconst < ijet->constituents().size(); iconst++) {
					theXConeGenJetConstPt     . push_back(ijet->constituents().at(iconst).pt());
					theXConeGenJetConstEta    . push_back(ijet->constituents().at(iconst).eta());
					theXConeGenJetConstPhi    . push_back(ijet->constituents().at(iconst).phi());
					theXConeGenJetConstEnergy . push_back(ijet->constituents().at(iconst).e());
			  }
			  /*if(DEBUG){
				for(unsigned int i=theXConeGenJetConstStartIndex.at(index); i <= theXConeGenJetConstEndIndex.at(index);i++){
					std::cout << "		const no: "<< i ;
					std::cout << "	pt: "<<theXConeGenJetConstPt.at(i);
					std::cout << "	eta: "<<theXConeGenJetConstEta.at(i);
					std::cout << "	phi: "<<theXConeGenJetConstPhi.at(i);
					std::cout << "	energy: "<<theXConeGenJetConstEnergy.at(i); 
					if(theXConeGenJetConstPt.at(i)<1e-50)std::cout << "	-----> GHOST!!"; 
					std::cout << std::endl;
				}
			  }*/
		  }
		  //collect constituent info - end

		}

		SetValue("theXConeGenJetPt",     theXConeGenJetPt);
		SetValue("theXConeGenJetEta",    theXConeGenJetEta);
		SetValue("theXConeGenJetPhi",    theXConeGenJetPhi);
		SetValue("theXConeGenJetEnergy", theXConeGenJetEnergy);
		SetValue("theXConeGenJetArea", theXConeGenJetArea);

		SetValue("theXConeGenJetConstStartIndex",     theXConeGenJetConstStartIndex);
		SetValue("theXConeGenJetConstEndIndex",     theXConeGenJetConstEndIndex);
		SetValue("theXConeGenJetConstPt",     theXConeGenJetConstPt);
		SetValue("theXConeGenJetConstEta",    theXConeGenJetConstEta);
		SetValue("theXConeGenJetConstPhi",    theXConeGenJetConstPhi);
		SetValue("theXConeGenJetConstEnergy", theXConeGenJetConstEnergy);

	}
	
	//Implementing XCone - end - Rizki

    return 0;
}

int XConeCalc::EndJob()
{
  return 0;
}
