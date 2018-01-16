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

//JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//JER 
#include "JetMETCorrections/Modules/interface/JetResolution.h"

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

  std::string MCL1JetPar;
  std::string MCL2JetPar;
  std::string MCL3JetPar;
  std::string JECunc_txtfile;

  std::string DataL1JetPar_BCD;
  std::string DataL2JetPar_BCD;
  std::string DataL3JetPar_BCD;
  std::string DataResJetPar_BCD;

  std::string DataL1JetPar_EF;
  std::string DataL2JetPar_EF;
  std::string DataL3JetPar_EF;
  std::string DataResJetPar_EF;

  std::string DataL1JetPar_G;
  std::string DataL2JetPar_G;
  std::string DataL3JetPar_G;
  std::string DataResJetPar_G;

  std::string DataL1JetPar_H;
  std::string DataL2JetPar_H;
  std::string DataL3JetPar_H;
  std::string DataResJetPar_H;

  std::string JER_txtfile;
  std::string JERSF_txtfile;
  
  //PUPPI JEC JEC
  std::vector<std::string>  jecPayloadsAK4chs;
  boost::shared_ptr<FactorizedJetCorrector>   JetCorrectorAK4chs;
  boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4chs;
  JME::JetResolution resolution_AK4chs;
  JME::JetResolutionScaleFactor resolution_SF_AK4chs;
  
  //PUPPI JEC JEC
  std::vector<std::string>  jecPayloadsAK4pup;
  boost::shared_ptr<FactorizedJetCorrector>   JetCorrectorAK4pup;
  boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4pup;  
  JME::JetResolution resolution_AK4pup;
  JME::JetResolutionScaleFactor resolution_SF_AK4pup;

  
  TRandom3 JERrand;
  
  bool JECup;
  bool JECdown;
  bool JERup;
  bool JERdown;

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
    std::cout << "XConeCalc: isMc = " << isMc << std::endl;

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

//    // Apply JEC?
//     if(mPset.exists("applyJEC")) applyJEC = mPset.getParameter<bool>("applyJEC");
//     else applyJEC = true; //default
//     std::cout << "XConeCalc: applyJEC = " << applyJEC << std::endl;

	//JECtxtfiles for MC
    if(mPset.exists("MCL1JetPar")) MCL1JetPar = mPset.getParameter<std::string>("MCL1JetPar");
    else MCL1JetPar = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L1FastJet";
    if(DEBUG) std::cout << "XConeCalc: Applying (MC) L1 Correction file: " << MCL1JetPar << std::endl;
    if(mPset.exists("MCL2JetPar")) MCL2JetPar = mPset.getParameter<std::string>("MCL2JetPar");
    else MCL2JetPar = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L2Relative";
    if(DEBUG) std::cout << "XConeCalc: Applying (MC) L2 Correction file: " << MCL2JetPar << std::endl;
    if(mPset.exists("MCL3JetPar")) MCL3JetPar = mPset.getParameter<std::string>("MCL3JetPar");
    else MCL3JetPar = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_L3Absolute";
    if(DEBUG) std::cout << "XConeCalc: Applying (MC) L3 Correction file: " << MCL3JetPar << std::endl;
    if(mPset.exists("JECunc_txtfile")) JECunc_txtfile = mPset.getParameter<std::string>("JECunc_txtfile");
    else JECunc_txtfile = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016V3_MC_Uncertainty";
    if(DEBUG) std::cout << "XConeCalc: Applying (MC) JECunc file: " << JECunc_txtfile << std::endl;

	//JER
    if(mPset.exists("JER_txtfile")) JER_txtfile = mPset.getParameter<std::string>("JER_txtfile");
    else JER_txtfile = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_PtResolution";
    if(DEBUG) std::cout << "XConeCalc: Applying (MC) JER file: " << JER_txtfile << std::endl;
    if(mPset.exists("JERSF_txtfile")) JERSF_txtfile = mPset.getParameter<std::string>("JERSF_txtfile");
    else JERSF_txtfile = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Spring16V10/Spring16_25nsV10_MC_SF";
    if(DEBUG) std::cout << "XConeCalc: Applying (MC) JERSF file: " << JERSF_txtfile << std::endl;	  


	//JEC and JER up/downs
    if(mPset.exists("JECup")) JECup = mPset.getParameter<bool>("JECup");
    else JECup = false;
    if(mPset.exists("JECdown")) JECdown = mPset.getParameter<bool>("JECdown");
    else JECdown = false;
    if(mPset.exists("JERup")) JERup = mPset.getParameter<bool>("JERup");
    else JERup = false;
    if(mPset.exists("JERdown")) JERdown = mPset.getParameter<bool>("JERdown");
    else JERdown = false;

	//JECpayloads/textfiles for DATA - BCD
    if(mPset.exists("DataL1JetPar_BCD")) DataL1JetPar_BCD = mPset.getParameter<std::string>("DataL1JetPar_BCD");
    else DataL1JetPar_BCD = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L1FastJet";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L1 Correction file: " << DataL1JetPar_BCD << std::endl;
    if(mPset.exists("DataL2JetPar_BCD")) DataL2JetPar_BCD = mPset.getParameter<std::string>("DataL2JetPar_BCD");
    else DataL2JetPar_BCD = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2Relative";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2 Correction file: " << DataL2JetPar_BCD << std::endl;
    if(mPset.exists("DataL3JetPar_BCD")) DataL3JetPar_BCD = mPset.getParameter<std::string>("DataL3JetPar_BCD");
    else DataL3JetPar_BCD = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L3Absolute";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L3 Correction file: " << DataL3JetPar_BCD << std::endl;
    if(mPset.exists("DataResJetPar_BCD")) DataResJetPar_BCD = mPset.getParameter<std::string>("DataResJetPar_BCD");
    else DataResJetPar_BCD = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016BCDV3_DATA_L2L3Residual";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2L3Res Correction file: " << DataResJetPar_BCD << std::endl;

	//JECpayloads/textfiles for DATA - EF
    if(mPset.exists("DataL1JetPar_EF")) DataL1JetPar_EF = mPset.getParameter<std::string>("DataL1JetPar_EF");
    else DataL1JetPar_EF = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016EFV3_DATA_L1FastJet";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L1 Correction file: " << DataL1JetPar_EF << std::endl;
    if(mPset.exists("DataL2JetPar_EF")) DataL2JetPar_EF = mPset.getParameter<std::string>("DataL2JetPar_EF");
    else DataL2JetPar_EF = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016EFV3_DATA_L2Relative";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2 Correction file: " << DataL2JetPar_EF << std::endl;
    if(mPset.exists("DataL3JetPar_EF")) DataL3JetPar_EF = mPset.getParameter<std::string>("DataL3JetPar_EF");
    else DataL3JetPar_EF = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016EFV3_DATA_L3Absolute";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L3 Correction file: " << DataL3JetPar_EF << std::endl;
    if(mPset.exists("DataResJetPar_EF")) DataResJetPar_EF = mPset.getParameter<std::string>("DataResJetPar_EF");
    else DataResJetPar_EF = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016EFV3_DATA_L2L3Residual";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2L3Res Correction file: " << DataResJetPar_EF << std::endl;

	//JECpayloads/textfiles for DATA - G
    if(mPset.exists("DataL1JetPar_G")) DataL1JetPar_G = mPset.getParameter<std::string>("DataL1JetPar_G");
    else DataL1JetPar_G = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016GV3_DATA_L1FastJet";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L1 Correction file: " << DataL1JetPar_G << std::endl;
    if(mPset.exists("DataL2JetPar_G")) DataL2JetPar_G = mPset.getParameter<std::string>("DataL2JetPar_G");
    else DataL2JetPar_G = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016GV3_DATA_L2Relative";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2 Correction file: " << DataL2JetPar_G << std::endl;
    if(mPset.exists("DataL3JetPar_G")) DataL3JetPar_G = mPset.getParameter<std::string>("DataL3JetPar_G");
    else DataL3JetPar_G = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016GV3_DATA_L3Absolute";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L3 Correction file: " << DataL3JetPar_G << std::endl;
    if(mPset.exists("DataResJetPar_G")) DataResJetPar_G = mPset.getParameter<std::string>("DataResJetPar_G");
    else DataResJetPar_G = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016GV3_DATA_L2L3Residual";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2L3Res Correction file: " << DataResJetPar_G << std::endl;

	//JECpayloads/textfiles for DATA - H
    if(mPset.exists("DataL1JetPar_H")) DataL1JetPar_H = mPset.getParameter<std::string>("DataL1JetPar_H");
    else DataL1JetPar_H = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016HV3_DATA_L1FastJet";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L1 Correction file: " << DataL1JetPar_H << std::endl;
    if(mPset.exists("DataL2JetPar_H")) DataL2JetPar_H = mPset.getParameter<std::string>("DataL2JetPar_H");
    else DataL2JetPar_H = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016HV3_DATA_L2Relative";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2 Correction file: " << DataL2JetPar_H << std::endl;
    if(mPset.exists("DataL3JetPar_H")) DataL3JetPar_H = mPset.getParameter<std::string>("DataL3JetPar_H");
    else DataL3JetPar_H = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016HV3_DATA_L3Absolute";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L3 Correction file: " << DataL3JetPar_H << std::endl;
    if(mPset.exists("DataResJetPar_H")) DataResJetPar_H = mPset.getParameter<std::string>("DataResJetPar_H");
    else DataResJetPar_H = "/uscms_data/d3/rsyarif/Fermilab2017/XConeinLJMet80x/CMSSW_8_0_25/src/LJMet/Com/data/Summer16RRV3/Summer16_23Sep2016HV3_DATA_L2L3Residual";
    if(DEBUG) std::cout << "XConeCalc: Applying (DATA) L2L3Res Correction file: " << DataResJetPar_H << std::endl;	
	
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

  	float		eventRho;

	//Implementing XCone - start - Rizki
	
	std::vector<fastjet::PseudoJet> FJConstituents; 
  	int		XConeNumJets_optimal=0;

	std::vector<double> RawXConeJetPt;
    std::vector<double> RawXConeJetEta;
    std::vector<double> RawXConeJetPhi;
    std::vector<double> RawXConeJetEnergy;
    std::vector<double> RawXConeJetArea;

	std::vector<double> theXConeJetPt;
    std::vector<double> theXConeJetEta;
    std::vector<double> theXConeJetPhi;
    std::vector<double> theXConeJetEnergy;
    std::vector<double> theXConeJetArea;
    std::vector<double> theXConeJetCorrScale;

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
    std::vector<double> theXConePUPPIJetCorrScale;
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
   if(isMc & doGenXCone){
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
	
	//SET UP Rho for JEC
    edm::Handle<double> rhoHandle;
    edm::InputTag rhoSrc_("fixedGridRhoFastjetAll", "");
    event.getByLabel(rhoSrc_, rhoHandle);
    double rho = std::max(*(rhoHandle.product()), 0.0);
	if (DEBUG) cout<<"SETUP: rho = "<<rho<<endl;	
	SetValue("eventRho",     rho);


	{ //make PFchs XCone jets - start

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
   
   // SET UP AK4chs JECpayloads (for XCone chs) 
   std::vector<JetCorrectorParameters> vParAK4chs;
   if(isMc){
	   jecPayloadsAK4chs.push_back(MCL1JetPar+"_AK4PFchs.txt");
	   jecPayloadsAK4chs.push_back(MCL2JetPar+"_AK4PFchs.txt");
	   jecPayloadsAK4chs.push_back(MCL3JetPar+"_AK4PFchs.txt");
	   jecPayloadsAK4chs.push_back(JECunc_txtfile+"_AK4PFchs.txt");
   }
   else{ //IMPLEMENT DATA JECpayloads for difference eras.
	   int iRun   = event.id().run();
	   if(iRun <= 276811){
	   		jecPayloadsAK4chs.push_back(DataL1JetPar_BCD+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL2JetPar_BCD+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL3JetPar_BCD+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataResJetPar_BCD+"_AK4PFchs.txt");
	   }
	   else if(iRun <= 278801){
	   		jecPayloadsAK4chs.push_back(DataL1JetPar_EF+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL2JetPar_EF+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL3JetPar_EF+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataResJetPar_EF+"_AK4PFchs.txt");	   
	   }
	   else if(iRun <= 280385){
	   		jecPayloadsAK4chs.push_back(DataL1JetPar_G+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL2JetPar_G+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL3JetPar_G+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataResJetPar_G+"_AK4PFchs.txt");	   
	   }
	   else{
	   		jecPayloadsAK4chs.push_back(DataL1JetPar_H+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL2JetPar_H+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataL3JetPar_H+"_AK4PFchs.txt");
	   		jecPayloadsAK4chs.push_back(DataResJetPar_H+"_AK4PFchs.txt");	   
	   }
   }
//    for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4chs.begin(), ipayloadEnd = jecPayloadsAK4chs.end(); ipayload != ipayloadEnd - 1; ++ipayload ) { //use this if want to skip last entry 
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4chs.begin(), ipayloadEnd = jecPayloadsAK4chs.end(); ipayload != ipayloadEnd; ++ipayload ) { //use this if use all entry
		
		if(isMc && (ipayload == ipayloadEnd - 1) ) continue; //for MC, skip JECuncer_txtfile at the last entry
   		if (DEBUG)cout<<"AK4chs JEC txt: "<<*ipayload<<endl;
   		JetCorrectorParameters pars(*ipayload);
   		vParAK4chs.push_back(pars);
   }

   JetCorrectorAK4chs   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4chs) );

   //Implementing JEC uncert for MC
   if(isMc) JetCorrUncertAK4chs  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty( jecPayloadsAK4chs.back() ) );

	//SET UP JER payload
    resolution_AK4chs = JME::JetResolution(JER_txtfile+"_AK4PFchs.txt");
    resolution_SF_AK4chs = JME::JetResolutionScaleFactor(JERSF_txtfile+"_AK4PFchs.txt");    


    if(DEBUG)std::cout << "---- " << N <<" XCone Jets ---- R = "<< XConeR << std::endl;
    int ConstituentStartIndex =0;
    for (std::vector<fastjet::PseudoJet>::const_iterator ijet = xcone_jetsA.begin(); ijet != xcone_jetsA.end(); ijet++) {
      int index = (int)(ijet-xcone_jetsA.begin());
      if(DEBUG)std::cout << "no. : " << index << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet constituents	: "<< ijet->constituents().size() << "	("<<ConstituentStartIndex <<" - "<< ConstituentStartIndex+ijet->constituents().size()-1 << ") "<<std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet pt			: "<< ijet->pt() << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet eta			: "<< ijet->eta() << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet phi			: "<< ijet->phi() << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet energy		: "<< ijet->e() << std::endl;
      if(DEBUG)std::cout << "      " <<  "Jet area			: "<< ijet->area() << std::endl;
      
      RawXConeJetPt     . push_back(ijet->pt());
      RawXConeJetEta    . push_back(ijet->eta());
      RawXConeJetPhi    . push_back(ijet->phi());
      RawXConeJetEnergy . push_back(ijet->e());
      RawXConeJetArea   . push_back(ijet->area());
      
      if(DEBUG) std::cout << "	--- Applying AK4CHS JEC correction to XCone" << std::endl;  
      //------------------------------------
      // Applying AK4CHS JEC correction to XCone
      //------------------------------------

      reco::Candidate::LorentzVector uncorrJet;
      uncorrJet.SetPx((double)ijet->px());
      uncorrJet.SetPy((double)ijet->py());
      uncorrJet.SetPz((double)ijet->pz());
      uncorrJet.SetE((double)ijet->e());
      if ( DEBUG ) cout << "   	-> before JEC pt,eta,phi,e	= " << uncorrJet.pt() << ",	" << uncorrJet.eta() << ",	" << uncorrJet.phi() << ",	" << uncorrJet.e() << endl;
      
      JetCorrectorAK4chs->setJetPt( ijet->pt() );
      JetCorrectorAK4chs->setJetEta ( ijet->eta() );
      JetCorrectorAK4chs->setJetE  ( ijet->e() );
      JetCorrectorAK4chs->setJetA  ( ijet->area() );
      JetCorrectorAK4chs->setRho   ( rho );
      double JECcorr = JetCorrectorAK4chs->getCorrection();

      reco::Candidate::LorentzVector corrJet = JECcorr * uncorrJet;
//       if ( DEBUG ) cout << "   -> after JEC pt,eta,phi,m = " << corrJet.pt() << ", " << corrJet.eta() << ", " << corrJet.phi() << ", " << corrJet.mass() << endl;
      if (DEBUG) cout << "		JECcorr : " << JECcorr << endl;
      if ( DEBUG ) cout << "   	-> after  JEC pt,eta,phi,e	= " << corrJet.pt() << ",	" << corrJet.eta() << ",	" << corrJet.phi() << ",	" << corrJet.e() << endl;
      
      double JERptscale = 1.0;
      double JECuncert = 1.0;
      if(isMc){

		  //------------------------------------
		  // Setting up AK4CHS JER for XCone (includes up/downs)
		  //------------------------------------
	  
		  Variation JERsystematic = Variation::NOMINAL;
		  if(JERup) JERsystematic = Variation::UP;
		  if(JERdown) JERsystematic = Variation::DOWN;

		  JME::JetParameters parameters;
		  parameters.setJetPt(corrJet.pt());
		  parameters.setJetEta(corrJet.eta());
		  parameters.setRho(rho);
		  double res = resolution_AK4chs.getResolution(parameters);
		  double factor = resolution_SF_AK4chs.getScaleFactor(parameters,JERsystematic) - 1;
	  
		  if (factor>0) {
			JERrand.SetSeed(abs(static_cast<int>(ijet->phi()*1e4)));
			JERptscale = max(0.0, JERrand.Gaus(corrJet.pt(),sqrt(factor*(factor+2))*res*corrJet.pt())/corrJet.pt());
		  }

		  //------------------------------------
		  // Setting up AK4CHS JEC uncertainty for XCone
		  //------------------------------------
	  
		  if ( JECup || JECdown ) {

				JetCorrUncertAK4chs->setJetEta(corrJet.eta());
				JetCorrUncertAK4chs->setJetPt(corrJet.pt()*JERptscale); //why does JECunc requires JER first? 

				if (JECup) { 
					try{
							JECuncert = JetCorrUncertAK4chs->getUncertainty(true);
					}
					catch(...){ // catch all exceptions. Jet Uncertainty tool throws when binning out of range
						std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
						std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
						std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
						JECuncert = 0.0;
				}
				JECuncert = 1 + JECuncert; 
		  }
				if (JECdown){ 
					try{
							JECuncert = JetCorrUncertAK4chs->getUncertainty(false);
					}
					catch(...){
						std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
						std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
						std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
						JECuncert = 0.0;
					}
						JECuncert = 1 - JECuncert; 
				}

				if (corrJet.pt()*JERptscale < 10.0 && JECup) JECuncert = 2.0;
				if (corrJet.pt()*JERptscale < 10.0 && JECdown) JECuncert = 0.01;

		  }

		  //------------------------------------
		  // Apply all jet corrections and uncertainties for XCone
		  //------------------------------------
	  
		  corrJet = JECuncert * JERptscale * JECcorr * uncorrJet;
		  if (DEBUG) cout << "		JERptscale : " << JERptscale << endl;
		  if (DEBUG) cout << "   	-> after JER + JEC + up/downs. pt,eta,phi,e	= " << corrJet.pt() << ",	" << corrJet.eta() << ",	" << corrJet.phi() << ",	" << corrJet.e() << endl;
		  if (DEBUG) cout << "		JECuncert : " << JECuncert << " (1 = Nominal) "<< endl;
	  
      }


      //------------------------------------
      // Save corrected XCone jets
      //------------------------------------

      theXConeJetPt     . push_back(corrJet.pt());
      theXConeJetEta    . push_back(corrJet.eta());
      theXConeJetPhi    . push_back(corrJet.phi());
      theXConeJetEnergy . push_back(corrJet.e());
      theXConeJetArea   . push_back(ijet->area());
      theXConeJetCorrScale   . push_back(JECuncert * JERptscale * JECcorr);

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

    SetValue("RawXConeJetPt",     RawXConeJetPt);
    SetValue("RawXConeJetEta",    RawXConeJetEta);
    SetValue("RawXConeJetPhi",    RawXConeJetPhi);
    SetValue("RawXConeJetEnergy", RawXConeJetEnergy);
    SetValue("RawXConeJetArea", RawXConeJetArea);

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
    
    } //make PFchs XCone jets - end

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
		
		// SET UP AK4pup JECpayloads (for XCone chs) 
		std::vector<JetCorrectorParameters> vParAK4pup;
		if(isMc){
		   jecPayloadsAK4pup.push_back(MCL1JetPar+"_AK4PFPuppi.txt");
		   jecPayloadsAK4pup.push_back(MCL2JetPar+"_AK4PFPuppi.txt");
		   jecPayloadsAK4pup.push_back(MCL3JetPar+"_AK4PFPuppi.txt");
		   jecPayloadsAK4pup.push_back(JECunc_txtfile+"_AK4PFPuppi.txt");
		}
		else{ //IMPLEMENT DATA JECpayloads for difference eras.
		   int iRun   = event.id().run();
		   if(iRun <= 276811){
				jecPayloadsAK4pup.push_back(DataL1JetPar_BCD+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL2JetPar_BCD+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL3JetPar_BCD+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataResJetPar_BCD+"_AK4PFPuppi.txt");
		   }
		   else if(iRun <= 278801){
				jecPayloadsAK4pup.push_back(DataL1JetPar_EF+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL2JetPar_EF+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL3JetPar_EF+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataResJetPar_EF+"_AK4PFPuppi.txt");	   
		   }
		   else if(iRun <= 280385){
				jecPayloadsAK4pup.push_back(DataL1JetPar_G+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL2JetPar_G+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL3JetPar_G+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataResJetPar_G+"_AK4PFPuppi.txt");	   
		   }
		   else{
				jecPayloadsAK4pup.push_back(DataL1JetPar_H+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL2JetPar_H+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataL3JetPar_H+"_AK4PFPuppi.txt");
				jecPayloadsAK4pup.push_back(DataResJetPar_H+"_AK4PFPuppi.txt");	   
		   }
		}
		
		for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4pup.begin(), ipayloadEnd = jecPayloadsAK4pup.end(); ipayload != ipayloadEnd; ++ipayload ) { //use this if use all entry
	
			if(isMc && (ipayload == ipayloadEnd - 1) ) continue; //for MC, skip JECuncer_txtfile at the last entry
			if (DEBUG)cout<<"AK4pup JEC txt: "<<*ipayload<<endl;
			JetCorrectorParameters pars(*ipayload);
			vParAK4pup.push_back(pars);
		}
		
		JetCorrectorAK4pup   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4pup) );
		
		//Implementing JEC uncert for MC
		if(isMc) JetCorrUncertAK4pup  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty( jecPayloadsAK4pup.back() ) );

		//SET UP JER payload
		resolution_AK4pup = JME::JetResolution(JER_txtfile+"_AK4PFPuppi.txt");
		resolution_SF_AK4pup = JME::JetResolutionScaleFactor(JERSF_txtfile+"_AK4PFPuppi.txt");    


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

		  if(DEBUG) std::cout << "	--- Applying AK4pup JEC correction to XCone" << std::endl;  
		  //------------------------------------
		  // Applying AK4pup JEC correction to XCone
		  //------------------------------------

		  reco::Candidate::LorentzVector uncorrJet;
		  uncorrJet.SetPx((double)ijet->px());
		  uncorrJet.SetPy((double)ijet->py());
		  uncorrJet.SetPz((double)ijet->pz());
		  uncorrJet.SetE((double)ijet->e());
		  if ( DEBUG ) cout << "   	-> before JEC pt,eta,phi,e	= " << uncorrJet.pt() << ",	" << uncorrJet.eta() << ",	" << uncorrJet.phi() << ",	" << uncorrJet.e() << endl;
	  
		  JetCorrectorAK4pup->setJetPt( ijet->pt() );
		  JetCorrectorAK4pup->setJetEta ( ijet->eta() );
		  JetCorrectorAK4pup->setJetE  ( ijet->e() );
		  JetCorrectorAK4pup->setJetA  ( ijet->area() );
		  JetCorrectorAK4pup->setRho   ( rho );
		  double JECcorr = JetCorrectorAK4pup->getCorrection();

		  reco::Candidate::LorentzVector corrJet = JECcorr * uncorrJet;
	//       if ( DEBUG ) cout << "   -> after JEC pt,eta,phi,m = " << corrJet.pt() << ", " << corrJet.eta() << ", " << corrJet.phi() << ", " << corrJet.mass() << endl;
		  if (DEBUG) cout << "		JECcorr : " << JECcorr << endl;
		  if ( DEBUG ) cout << "   	-> after  JEC pt,eta,phi,e	= " << corrJet.pt() << ",	" << corrJet.eta() << ",	" << corrJet.phi() << ",	" << corrJet.e() << endl;
		  
		  double JERptscale = 1.0;
		  double JECuncert = 1.0;
		  if(isMc){

			  //------------------------------------
			  // Setting up AK4pup JER for XCone (includes up/downs)
			  //------------------------------------
	  
			  Variation JERsystematic = Variation::NOMINAL;
			  if(JERup) JERsystematic = Variation::UP;
			  if(JERdown) JERsystematic = Variation::DOWN;

			  JME::JetParameters parameters;
			  parameters.setJetPt(corrJet.pt());
			  parameters.setJetEta(corrJet.eta());
			  parameters.setRho(rho);
			  double res = resolution_AK4pup.getResolution(parameters);
			  double factor = resolution_SF_AK4pup.getScaleFactor(parameters,JERsystematic) - 1;
	  
			  if (factor>0) {
				JERrand.SetSeed(abs(static_cast<int>(ijet->phi()*1e4)));
				JERptscale = max(0.0, JERrand.Gaus(corrJet.pt(),sqrt(factor*(factor+2))*res*corrJet.pt())/corrJet.pt());
			  }

			  //------------------------------------
			  // Setting up AK4pup JEC uncertainty for XCone
			  //------------------------------------
			  	  
			  if ( JECup || JECdown ) {

					JetCorrUncertAK4pup->setJetEta(corrJet.eta());
					JetCorrUncertAK4pup->setJetPt(corrJet.pt()*JERptscale); //why does JECunc requires JER first? 

					if (JECup) { 
						try{
								JECuncert = JetCorrUncertAK4pup->getUncertainty(true);
						}
						catch(...){ // catch all exceptions. Jet Uncertainty tool throws when binning out of range
							std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
							std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
							std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
							JECuncert = 0.0;
					}
					JECuncert = 1 + JECuncert; 
			  }
					if (JECdown){ 
						try{
								JECuncert = JetCorrUncertAK4pup->getUncertainty(false);
						}
						catch(...){
							std::cout << mLegend << "WARNING! Exception thrown by JetCorrectionUncertainty!" << std::endl;
							std::cout << mLegend << "WARNING! Possibly, trying to correct a jet/MET outside correction range." << std::endl;
							std::cout << mLegend << "WARNING! Jet/MET will remain uncorrected." << std::endl;
							JECuncert = 0.0;
						}
							JECuncert = 1 - JECuncert; 
					}

					if (corrJet.pt()*JERptscale < 10.0 && JECup) JECuncert = 2.0;
					if (corrJet.pt()*JERptscale < 10.0 && JECdown) JECuncert = 0.01;

			  }

			  //------------------------------------
			  // Apply all jet corrections and uncertainties for XCone
			  //------------------------------------
	  
			  corrJet = JECuncert * JERptscale * JECcorr * uncorrJet;
			  if (DEBUG) cout << "		JERptscale : " << JERptscale << endl;
			  if (DEBUG) cout << "   	-> after JER + JEC + up/downs. pt,eta,phi,e	= " << corrJet.pt() << ",	" << corrJet.eta() << ",	" << corrJet.phi() << ",	" << corrJet.e() << endl;
			  if (DEBUG) cout << "		JECuncert : " << JECuncert << " (1 = Nominal) "<< endl;
	  
		  }


		  //------------------------------------
		  // Save corrected XCone jets
		  //------------------------------------

  
		  theXConePUPPIJetPt     . push_back(ijet->pt());
		  theXConePUPPIJetEta    . push_back(ijet->eta());
		  theXConePUPPIJetPhi    . push_back(ijet->phi());
		  theXConePUPPIJetEnergy . push_back(ijet->e());
		  theXConePUPPIJetArea   . push_back(ijet->area());
		  theXConePUPPIJetCorrScale   . push_back(JECuncert * JERptscale * JECcorr);

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

	//Clear JEC vector str
	jecPayloadsAK4chs  .clear();	
	jecPayloadsAK4pup  .clear();	

    return 0;
}

int XConeCalc::EndJob()
{
  return 0;
}
