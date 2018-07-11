
/*
 Calculator for substructure variables
 
 Author: Joshua Swanson, 2014. Updated Julie Hogan, 2015. 
 */

#include <iostream>
#include <limits>   // std::numeric_limits
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TF1.h"

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/BaseEventSelector.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "LJMet/Com/interface/HTTTopJetTagInfo.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

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
  edm::InputTag genParticles_it;
  edm::InputTag httTagInfo_it;
  std::string bDiscriminant;
  std::string tagInfo;
  double kappa;
  bool killHF;
  bool doNewJEC;
  bool isMc;
  bool JECup;
  bool JECdn;
  bool JERup;
  bool JERdn;
  bool useL2L3Mass;
  std::string puppiCorrPath;
  std::string MCL3;
  std::string MCL2;
  std::string MCSF;
  std::string DataL3;
  std::string DataL2;
  std::string DataL2L3;
  JME::JetResolutionScaleFactor resolution_SF;
  FactorizedJetCorrector *jecak8;
  JetCorrectorParameters *L3JetParAK8MC;
  JetCorrectorParameters *L2JetParAK8MC;
  JetCorrectorParameters *ResJetParAK8MC; 
  JetCorrectorParameters *L3JetParAK8BCD;
  JetCorrectorParameters *L2JetParAK8BCD;
  JetCorrectorParameters *ResJetParAK8BCD; 
  JetCorrectorParameters *L3JetParAK8EF;
  JetCorrectorParameters *L2JetParAK8EF;
  JetCorrectorParameters *ResJetParAK8EF; 
  JetCorrectorParameters *L3JetParAK8G;
  JetCorrectorParameters *L2JetParAK8G;
  JetCorrectorParameters *ResJetParAK8G; 
  JetCorrectorParameters *L3JetParAK8H;
  JetCorrectorParameters *L2JetParAK8H;
  JetCorrectorParameters *ResJetParAK8H; 
  JetCorrectionUncertainty *jecUnc;
  std::vector<JetCorrectorParameters> vParAK8MC;
  std::vector<JetCorrectorParameters> vParAK8BCD;
  std::vector<JetCorrectorParameters> vParAK8EF;
  std::vector<JetCorrectorParameters> vParAK8G;
  std::vector<JetCorrectorParameters> vParAK8H;
  TRandom3 JERrand;
  TF1 *puppisd_corrGEN;
  TF1 *puppisd_corrRECO_cen;
  TF1 *puppisd_corrRECO_for;
};

static int reg = LjmetFactory::GetInstance()->Register(new JetSubCalc(), "JetSubCalc");

struct sortclass {
  //  bool operator() (double i,double j) { return (i>j);}
  bool operator() (pat::Jet i, pat::Jet j){ return (i.pt() > j.pt());}
} ptsorter;


JetSubCalc::JetSubCalc()
{
}

JetSubCalc::~JetSubCalc()
{
}

int JetSubCalc::BeginJob()
{
    if (mPset.exists("slimmedJetColl")){ slimmedJetColl_it = mPset.getParameter<edm::InputTag>("slimmedJetColl");}
    else slimmedJetColl_it = edm::InputTag("slimmedJets");
    
    if (mPset.exists("slimmedJetsAK8Coll")) slimmedJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("slimmedJetsAK8Coll");
    else slimmedJetsAK8Coll_it = edm::InputTag("slimmedJetsAK8");

    if (mPset.exists("genParticles")) genParticles_it = mPset.getParameter<edm::InputTag>("genParticles");
    else                              genParticles_it = edm::InputTag("prunedGenParticles");

    if (mPset.exists("bDiscriminant")) bDiscriminant = mPset.getParameter<std::string>("bDiscriminant");
    else bDiscriminant = "pfDeepCSVJetTags:probb";
    
    std::cout << " JetSubCalc Bdisc = " << bDiscriminant << std::endl;
    
    if (mPset.exists("tagInfo")) tagInfo = mPset.getParameter<std::string>("tagInfo");
    else tagInfo = "CATop";

    if(mPset.exists("kappa")) kappa = mPset.getParameter<double>("kappa");
    else kappa = 0.5;

    if(mPset.exists("killHF")) killHF = mPset.getParameter<bool>("killHF");
    else killHF = false;

    if(mPset.exists("doNewJEC")) doNewJEC = mPset.getParameter<bool>("doNewJEC");
    else doNewJEC = false;

    cout << "JetSubCalc: doNewJEC = " << doNewJEC << ", killHF = " << killHF << endl;

    if(mPset.exists("useL2L3Mass")) useL2L3Mass = mPset.getParameter<bool>("useL2L3Mass");
    else useL2L3Mass = false;

    if(mPset.exists("puppiCorrPath")) puppiCorrPath = mPset.getParameter<std::string>("puppiCorrPath");
    else puppiCorrPath = "CMSSW_8_0_26_patch1/src/LJMet/Com/PuppiSoftdropMassCorr/weights/puppiCorr.root";

    if(mPset.exists("isMc")) isMc = mPset.getParameter<bool>("isMc");
    else isMc = false;

    if(mPset.exists("JECup")) JECup = mPset.getParameter<bool>("JECup");
    else JECup = false;

    if(mPset.exists("JECdown")) JECdn = mPset.getParameter<bool>("JECdown");
    else JECdn = false;

    if(mPset.exists("JERup")) JERup = mPset.getParameter<bool>("JERup");
    else JERup = false;

    if(mPset.exists("JERdown")) JERdn = mPset.getParameter<bool>("JERdown");
    else JERdn = false;

    if(useL2L3Mass){
      cout << "JetSubCalc: using L2+L3 corrected groomed masses" << endl;
      if(mPset.exists("MCL2JetParAK8")) MCL2 = mPset.getParameter<std::string>("MCL2JetParAK8");
      else MCL2 = "/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_MC_L1FastJet";
      if(mPset.exists("MCL3JetParAK8")) MCL3 = mPset.getParameter<std::string>("MCL3JetParAK8");
      else MCL3 = "/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_MC_L1FastJet";
      if(mPset.exists("MCSF")) MCSF = mPset.getParameter<std::string>("MCSF");
      else MCSF = "/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_MC_L1FastJet";

      if(mPset.exists("DataL2JetParAK8")) DataL2 = mPset.getParameter<std::string>("DataL2JetParAK8");
      else DataL2 = "/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_Data_L1FastJet";
      if(mPset.exists("DataL3JetParAK8")) DataL3 = mPset.getParameter<std::string>("DataL3JetParAK8");
      else DataL3 = "/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_Data_L1FastJet";
      if(mPset.exists("DataL2L3JetParAK8")) DataL2L3 = mPset.getParameter<std::string>("DataL2L3JetParAK8");
      else DataL2L3 = "/uscms_data/d3/jmanagan/CMSSW_7_4_14/src/LJMet/Com/data/Summer15_25nsV5_Data_L1FastJet";
 
      std::string strL2BCD = DataL2;
      std::string strL2EF = strL2BCD; boost::replace_first(strL2EF,"BCD","EF");
      std::string strL2G = strL2BCD; boost::replace_first(strL2G,"BCD","G");
      std::string strL2H = strL2BCD; boost::replace_first(strL2H,"BCD","H");  

      std::string strL3BCD = DataL3;
      std::string strL3EF = strL3BCD; boost::replace_first(strL3EF,"BCD","EF");
      std::string strL3G = strL3BCD; boost::replace_first(strL3G,"BCD","G");
      std::string strL3H = strL3BCD; boost::replace_first(strL3H,"BCD","H");  

      std::string strL2L3BCD = DataL2L3;
      std::string strL2L3EF = strL2L3BCD; boost::replace_first(strL2L3EF,"BCD","EF");
      std::string strL2L3G = strL2L3BCD; boost::replace_first(strL2L3G,"BCD","G");
      std::string strL2L3H = strL2L3BCD; boost::replace_first(strL2L3H,"BCD","H");  
       
      if(isMc){
	L3JetParAK8MC  = new JetCorrectorParameters(MCL3);
	L2JetParAK8MC  = new JetCorrectorParameters(MCL2);
	vParAK8MC.push_back(*L2JetParAK8MC);
	vParAK8MC.push_back(*L3JetParAK8MC);
      }else{
	ResJetParAK8BCD = new JetCorrectorParameters(DataL2L3); 
	L3JetParAK8BCD  = new JetCorrectorParameters(DataL3);
	L2JetParAK8BCD  = new JetCorrectorParameters(DataL2);
	vParAK8BCD.push_back(*L2JetParAK8BCD);
	vParAK8BCD.push_back(*L3JetParAK8BCD);
	vParAK8BCD.push_back(*ResJetParAK8BCD);

	ResJetParAK8EF = new JetCorrectorParameters(strL2L3EF); 
	L3JetParAK8EF  = new JetCorrectorParameters(strL3EF);
	L2JetParAK8EF  = new JetCorrectorParameters(strL2EF);
	vParAK8EF.push_back(*L2JetParAK8EF);
	vParAK8EF.push_back(*L3JetParAK8EF);
	vParAK8EF.push_back(*ResJetParAK8EF);

	ResJetParAK8G = new JetCorrectorParameters(strL2L3G); 
	L3JetParAK8G  = new JetCorrectorParameters(strL3G);
	L2JetParAK8G  = new JetCorrectorParameters(strL2G);
	vParAK8G.push_back(*L2JetParAK8G);
	vParAK8G.push_back(*L3JetParAK8G);
	vParAK8G.push_back(*ResJetParAK8G);

	ResJetParAK8H = new JetCorrectorParameters(strL2L3H); 
	L3JetParAK8H  = new JetCorrectorParameters(strL3H);
	L2JetParAK8H  = new JetCorrectorParameters(strL2H);
	vParAK8H.push_back(*L2JetParAK8H);
	vParAK8H.push_back(*L3JetParAK8H);
	vParAK8H.push_back(*ResJetParAK8H);
      }
      
      if(isMc) jecak8 = new FactorizedJetCorrector(vParAK8MC);
      jecUnc = new JetCorrectionUncertainty(mPset.getParameter<std::string>("UncertaintyAK8"));

      resolution_SF = JME::JetResolutionScaleFactor(MCSF);    
    }

    TFile* file = TFile::Open(puppiCorrPath.c_str(),"READ");
    puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
    puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
    puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");
    file->Close();
    

    return 0;
}

int JetSubCalc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector)
{
    // ----- Get AK4 jet objects from the selector -----
    // This is updated -- original version used all AK4 jets without selection 
    std::vector<pat::Jet>                       const & theJets = selector->GetSelectedCorrJets();
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
       
    // Discriminator for the MVA PileUp id.
    // NOTE: Training used is for ak5PFJetsCHS in CMSSW 5.3.X and Run 1 pileup
    std::vector<double> theJetPileupJetId;
    
    //Identity
    std::vector<int> theJetIndex;
    std::vector<int> theJetnDaughters;
    
    //Daughter four std::vector and index
    std::vector<double> theJetDaughterPt;
    std::vector<double> theJetDaughterEta;
    std::vector<double> theJetDaughterPhi;
    std::vector<double> theJetDaughterEnergy;
    
    std::vector <double> theJetCSVb;
    std::vector <double> theJetCSVbb;
    std::vector <double> theJetCSVc;
    std::vector <double> theJetCSVudsg;

    std::vector<int> theJetDaughterMotherIndex;
    std::vector<int> theJetCSVLSubJets;
    std::vector<int> theJetCSVMSubJets;
    std::vector<int> theJetCSVTSubJets;

    std::vector<int> theJetPFlav;
    std::vector<int> theJetHFlav;
    std::vector<int> theJetBTag;
    std::vector<int> theJetBTag_bSFup;
    std::vector<int> theJetBTag_bSFdn;
    std::vector<int> theJetBTag_lSFup;
    std::vector<int> theJetBTag_lSFdn;

    std::vector<int> maxProb;
    
    double thePileupJetId;

    for (std::vector<pat::Jet>::const_iterator ijet = theJets.begin(); ijet != theJets.end(); ijet++) {
      int index = (int)(ijet-theJets.begin());

      if(killHF && fabs(ijet->eta()) > 2.4) continue;
      
      thePileupJetId = -std::numeric_limits<double>::max();
      thePileupJetId = (double)ijet->userFloat("pileupJetId:fullDiscriminant");
      theJetPileupJetId.push_back(thePileupJetId);

      theJetPt     . push_back(ijet->pt());
      theJetEta    . push_back(ijet->eta());
      theJetPhi    . push_back(ijet->phi());
      theJetEnergy . push_back(ijet->energy());
      
      theJetCSVb.    push_back(ijet->bDiscriminator( bDiscriminant ));
      theJetCSVbb.   push_back(ijet->bDiscriminator( "pfDeepCSVJetTags:probbb" ));
      theJetCSVc.    push_back(ijet->bDiscriminator( "pfDeepCSVJetTags:probc" ));
      theJetCSVudsg. push_back(ijet->bDiscriminator( "pfDeepCSVJetTags:probudsg" ));

      theJetBTag.push_back(selector->isJetTagged(*ijet, event, true, 0));
      theJetBTag_bSFup.push_back(selector->isJetTagged(*ijet, event, true, 1));
      theJetBTag_bSFdn.push_back(selector->isJetTagged(*ijet, event, true, 2));
      theJetBTag_lSFup.push_back(selector->isJetTagged(*ijet, event, true, 3));
      theJetBTag_lSFdn.push_back(selector->isJetTagged(*ijet, event, true, 4));
      theJetPFlav.push_back(abs(ijet->partonFlavour()));
      theJetHFlav.push_back(abs(ijet->hadronFlavour()));

      theJetIndex.push_back(index);
      theJetnDaughters.push_back((int)ijet->numberOfDaughters());
      
      //HT
      theJetHT += ijet->pt(); 

    }
    
    double leading_pt = -999.;
    for (std::vector<pat::Jet>::const_iterator ijet = theJets.begin(); ijet != theJets.end(); ijet++) { 

      if(killHF && fabs(ijet->eta()) > 2.4) continue;
      if (ijet->pt() > leading_pt){leading_pt = ijet->pt();}
    }

    double second_leading_pt =-999.;
    for (std::vector<pat::Jet>::const_iterator ijet = theJets.begin(); ijet != theJets.end(); ijet++) { 

      if(killHF && fabs(ijet->eta()) > 2.4) continue;
      if(ijet->pt() > second_leading_pt && ijet->pt() < leading_pt){
	second_leading_pt = ijet->pt();
      }
    }

    SetValue("theJetPt",     theJetPt);
    SetValue("theJetEta",    theJetEta);
    SetValue("theJetPhi",    theJetPhi);
    SetValue("theJetEnergy", theJetEnergy);

    SetValue("theJetCSVb",    theJetCSVb);
    SetValue("theJetCSVbb",   theJetCSVbb);
    SetValue("theJetCSVc",    theJetCSVc);
    SetValue("theJetCSVudsg", theJetCSVudsg);

    SetValue("theJetPFlav",  theJetPFlav);
    SetValue("theJetHFlav",  theJetHFlav);
    SetValue("theJetBTag",   theJetBTag);
    SetValue("theJetBTag_bSFup",   theJetBTag_bSFup);
    SetValue("theJetBTag_bSFdn",   theJetBTag_bSFdn);
    SetValue("theJetBTag_lSFup",   theJetBTag_lSFup);
    SetValue("theJetBTag_lSFdn",   theJetBTag_lSFdn);

    SetValue("theJetHT", theJetHT);    
    SetValue("theJetLeadPt", leading_pt);
    SetValue("theJetSubLeadPt", second_leading_pt);

    SetValue("theJetPileupJetId", theJetPileupJetId);
    SetValue("theJetnDaughters", theJetnDaughters);

    // Load in AK8 jets (no selection performed on these)
    //edm::Handle<std::vector<pat::Jet> > theAK8Jets;
    //event.getByLabel(slimmedJetsAK8Coll_it, theAK8Jets);
    std::vector<pat::Jet> const & theAK8Jets = selector->GetSelectedCorrJets_AK8();

    // Four std::vector
    std::vector<double> theJetAK8Pt;
    std::vector<double> theJetAK8Eta;
    std::vector<double> theJetAK8Phi;
    std::vector<double> theJetAK8Energy;
    std::vector<double> theJetAK8CSV;
    std::vector<double> theJetAK8DoubleB;
    std::vector<double> theJetAK8JetCharge;
    std::vector<double> theJetAK8GenPt;
    std::vector<double> theJetAK8GenDR;
    std::vector<double> theJetAK8GenMass;

    // CHS values
    std::vector<double> theJetAK8CHSPt;
    std::vector<double> theJetAK8CHSEta;
    std::vector<double> theJetAK8CHSPhi;
    std::vector<double> theJetAK8CHSMass;
    std::vector<double> theJetAK8CHSTau1;
    std::vector<double> theJetAK8CHSTau2;
    std::vector<double> theJetAK8CHSTau3;

    std::vector<double> theJetAK8SoftDropRaw;    
    std::vector<double> theJetAK8SoftDropCorr;    
    std::vector<double> theJetAK8CHSSoftDropMass;
    std::vector<double> theJetAK8CHSPrunedMass;

    std::vector<double> theJetAK8SoftDrop_JMSup;    
    std::vector<double> theJetAK8SoftDrop_JMSdn;    
    std::vector<double> theJetAK8SoftDrop_JMRup;
    std::vector<double> theJetAK8SoftDrop_JMRdn;    

    std::vector<int>    theJetAK8SJIndex;
    std::vector<int>    theJetAK8SJSize;

    
    // Pruned, trimmed and filtered masses available
    std::vector<double> theJetAK8SoftDrop;

    // n-subjettiness variables tau1, tau2, and tau3 available
    std::vector<double> theJetAK8NjettinessTau1;
    std::vector<double> theJetAK8NjettinessTau2;
    std::vector<double> theJetAK8NjettinessTau3;

    std::vector<double> theJetAK8Mass;
    std::vector<int>    theJetAK8Index;
    std::vector<int>    theJetAK8nDaughters;
    
    std::vector<double> theJetAK8SDSubjetPt;
    std::vector<double> theJetAK8SDSubjetEta;
    std::vector<double> theJetAK8SDSubjetPhi;
    std::vector<double> theJetAK8SDSubjetMass;
    std::vector<double> theJetAK8SDSubjetCSVb;
    std::vector<double> theJetAK8SDSubjetCSVc;
    std::vector<double> theJetAK8SDSubjetCSVudsg;
    std::vector<double> theJetAK8SDSubjetCSVbb;
    std::vector<int>    theJetAK8SDSubjetHFlav;
    std::vector<double> theJetAK8SDSubjetBTag;
    std::vector<double> theJetAK8SDSubjetDR;
    std::vector<int> theJetAK8SDSubjetIndex;
    std::vector<int> theJetAK8SDSubjetSize;
    std::vector<int> theJetAK8SDSubjetNDeepCSVL;
    std::vector<int> theJetAK8SDSubjetNDeepCSVMSF;
    std::vector<int> theJetAK8SDSubjetNDeepCSVM_lSFup;
    std::vector<int> theJetAK8SDSubjetNDeepCSVM_lSFdn;
    std::vector<int> theJetAK8SDSubjetNDeepCSVM_bSFup;
    std::vector<int> theJetAK8SDSubjetNDeepCSVM_bSFdn;

    std::vector<double> theJetAK8SoftDropn2b1;
    std::vector<double> theJetAK8SoftDropn3b1;
    std::vector<double> theJetAK8SoftDropn2b2;
    std::vector<double> theJetAK8SoftDropn3b2;
    
    double topMass, minMass, jetCharge;
    int nSubJets;
    double theSoftDrop;
    double theCHSPrunedMass,theCHSSoftDropMass;
    double theNjettinessTau1, theNjettinessTau2, theNjettinessTau3;
    double theCHSTau1, theCHSTau2, theCHSTau3;
    double theSoftDropn2b1, theSoftDropn3b1, theSoftDropn2b2, theSoftDropn3b2;
   
    double SDsubjetPt;
    double SDsubjetEta; 
    double SDsubjetPhi;      
    double SDsubjetMass;      
    double SDsubjetDeepCSVb;
    double SDsubjetDeepCSVc;
    double SDsubjetDeepCSVudsg;
    double SDsubjetDeepCSVbb;
    int    SDsubjetHFlav;
    double SDsubjetBTag;     
    double SDdeltaRsubjetJet; 

    int nSDSubJets;
    int nSDSubsDeepCSVL;
    int nSDSubsDeepCSVM;
    int nSDSubsDeepCSVMSF;
    int nSDSubsDeepCSVM_bSFup;
    int nSDSubsDeepCSVM_bSFdn;
    int nSDSubsDeepCSVM_lSFup;
    int nSDSubsDeepCSVM_lSFdn;
    int SDSubJetIndex;

    double CHSsubjetPt;
    double CHSsubjetEta; 
    double CHSsubjetPhi;      
    double CHSsubjetMass;      
    double CHSsubjetBdisc;     
    int    CHSsubjetHFlav;
    double CHSsubjetBTag;     

    int nCHSSubJets;
    int nCHSSubsCSVL;
    int nCHSSubsCSVM;
    int nCHSSubsCSVMSF;
    int nCHSSubsCSVM_bSFup;
    int nCHSSubsCSVM_bSFdn;
    int nCHSSubsCSVM_lSFup;
    int nCHSSubsCSVM_lSFdn;
    int CHSSubJetIndex;

    int iRun   = event.id().run();
    if(!isMc){
      delete jecak8;
      if(iRun <= 276811){
	jecak8 = new FactorizedJetCorrector(vParAK8BCD);
      }
      else if(iRun <= 278801){ 
	jecak8 = new FactorizedJetCorrector(vParAK8EF);
      }
      else if(iRun <= 280385){
	jecak8 = new FactorizedJetCorrector(vParAK8G);
      }
      else{
	jecak8 = new FactorizedJetCorrector(vParAK8H); 
      }
    }

    //for (std::vector<pat::Jet>::const_iterator ijet = theAK8Jets->begin(); ijet != theAK8Jets->end(); ijet++) {
    for (std::vector<pat::Jet>::const_iterator ii = theAK8Jets.begin(); ii != theAK8Jets.end(); ii++){
      int index = (int)(ii-theAK8Jets.begin());

      if (ii->pt() < 200) continue;

      pat::Jet corrak8;
      corrak8 = *ii;

      if(killHF && fabs(corrak8.eta()) > 2.4) continue;

      theJetAK8Pt    .push_back(corrak8.pt());
      theJetAK8Eta   .push_back(corrak8.eta());
      theJetAK8Phi   .push_back(corrak8.phi());
      theJetAK8Energy.push_back(corrak8.energy());
      theJetAK8Mass  .push_back(corrak8.mass());

      double theCHSPt = -std::numeric_limits<double>::max();
      double theCHSEta = -std::numeric_limits<double>::max();
      double theCHSPhi = -std::numeric_limits<double>::max();
      double theCHSMass = -std::numeric_limits<double>::max();
      theCHSPt = corrak8.userFloat("ak8PFJetsCHSValueMap:pt");  
      theCHSEta = corrak8.userFloat("ak8PFJetsCHSValueMap:eta"); 
      theCHSPhi = corrak8.userFloat("ak8PFJetsCHSValueMap:phi"); 
      theCHSMass = corrak8.userFloat("ak8PFJetsCHSValueMap:mass");
      theJetAK8CHSPt.push_back(theCHSPt);
      theJetAK8CHSEta.push_back(theCHSEta);
      theJetAK8CHSPhi.push_back(theCHSPhi);
      theJetAK8CHSMass.push_back(theCHSMass);

      double genDR = -99;
      double genpt = -99;
      double genmass = -99;
      TLorentzVector ak8jet;
      ak8jet.SetPtEtaPhiE(corrak8.pt(),corrak8.eta(),corrak8.phi(),corrak8.energy());
      const reco::GenJet * genJet = corrak8.genJet();
      if(genJet){
	TLorentzVector genP4;
	genP4.SetPtEtaPhiE(genJet->pt(),genJet->eta(),genJet->phi(),genJet->energy());
	genDR = ak8jet.DeltaR(genP4);	
	genpt = genJet->pt();
	genmass = genJet->mass();
      }
      theJetAK8GenPt.push_back(genpt);
      theJetAK8GenMass.push_back(genmass);
      theJetAK8GenDR.push_back(genDR);

      theCHSPrunedMass   = -std::numeric_limits<double>::max();
      theCHSSoftDropMass = -std::numeric_limits<double>::max();
      theSoftDrop = -std::numeric_limits<double>::max();
	
      theCHSPrunedMass   = (double)corrak8.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass");
      theCHSSoftDropMass = (double)corrak8.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass");
      theSoftDrop = (double)corrak8.userFloat("ak8PFJetsPuppiSoftDropMass");		       	

      theNjettinessTau1 = std::numeric_limits<double>::max();
      theNjettinessTau2 = std::numeric_limits<double>::max();
      theNjettinessTau3 = std::numeric_limits<double>::max();
      theNjettinessTau1 = (double)corrak8.userFloat("NjettinessAK8Puppi:tau1");
      theNjettinessTau2 = (double)corrak8.userFloat("NjettinessAK8Puppi:tau2");
      theNjettinessTau3 = (double)corrak8.userFloat("NjettinessAK8Puppi:tau3");

      theCHSTau1 = std::numeric_limits<double>::max();
      theCHSTau2 = std::numeric_limits<double>::max();
      theCHSTau3 = std::numeric_limits<double>::max();
      theCHSTau1 = (double)corrak8.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1");
      theCHSTau2 = (double)corrak8.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2");
      theCHSTau3 = (double)corrak8.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3");

      theSoftDropn2b1 = std::numeric_limits<double>::max();
      theSoftDropn3b1 = std::numeric_limits<double>::max();
      theSoftDropn2b2 = std::numeric_limits<double>::max();
      theSoftDropn3b2 = std::numeric_limits<double>::max();
      theSoftDropn2b1 = (double)corrak8.userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2");
      theSoftDropn3b1 = (double)corrak8.userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3");
      theSoftDropn2b2 = (double)corrak8.userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2");
      theSoftDropn3b2 = (double)corrak8.userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN3");

      theJetAK8CSV.push_back(corrak8.bDiscriminator( bDiscriminant ));
      theJetAK8DoubleB.push_back(corrak8.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
      
      theJetAK8CHSPrunedMass.push_back(theCHSPrunedMass); // JEC only
      theJetAK8CHSSoftDropMass.push_back(theCHSSoftDropMass); // JEC only
      theJetAK8SoftDrop.push_back(theSoftDrop);
      
      theJetAK8CHSTau1.push_back(theCHSTau1);
      theJetAK8CHSTau2.push_back(theCHSTau2);
      theJetAK8CHSTau3.push_back(theCHSTau3);
      theJetAK8NjettinessTau1.push_back(theNjettinessTau1);
      theJetAK8NjettinessTau2.push_back(theNjettinessTau2);
      theJetAK8NjettinessTau3.push_back(theNjettinessTau3);

      theJetAK8SoftDropn2b1.push_back(theSoftDropn2b1);
      theJetAK8SoftDropn2b2.push_back(theSoftDropn2b2);
      theJetAK8SoftDropn3b1.push_back(theSoftDropn3b1);
      theJetAK8SoftDropn3b2.push_back(theSoftDropn3b2);

      theJetAK8nDaughters.push_back((int)corrak8.numberOfDaughters());
      theJetAK8Index.push_back(index);

      //JetCharge calculation

      reco::Jet::Constituents constituents = corrak8.getJetConstituents();

      double sumWeightedCharge = 0.0;
      int con_charge = 0;
      double con_pt = 0.0;
      for(auto constituentItr=constituents.begin(); constituentItr!=constituents.end(); ++constituentItr){
	edm::Ptr<reco::Candidate> constituent=*constituentItr;

	con_charge = (int)constituent->charge();
	con_pt     = (double)constituent->pt();
	
	sumWeightedCharge = sumWeightedCharge + ( con_charge * pow(con_pt,kappa) );
	
      }

      jetCharge  = 1.0/( pow( (corrak8.pt()), kappa) ) * sumWeightedCharge;

      theJetAK8JetCharge.push_back(jetCharge);

      // Get Soft drop subjets for subjet b-tagging
      SDSubJetIndex = (int)theJetAK8SDSubjetPt.size();
      nSDSubJets  = std::numeric_limits<int>::min();
      nSDSubsDeepCSVL = 0;
      nSDSubsDeepCSVMSF = 0;
      nSDSubsDeepCSVM_bSFup = 0;
      nSDSubsDeepCSVM_bSFdn = 0;
      nSDSubsDeepCSVM_lSFup = 0;
      nSDSubsDeepCSVM_lSFdn = 0;
      double maxSubCSV = 0;

      TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
      auto const & sdSubjets = corrak8.subjets("SoftDropPuppi");
      nSDSubJets = (int)sdSubjets.size();
      for ( auto const & it : sdSubjets ) {

	puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
	puppi_softdrop+=puppi_softdrop_subjet;

	pat::Jet corrsubjet;
	if(doNewJEC){
	  corrsubjet = selector->correctJetReturnPatJet(*it, event, false);
	}else{
	  corrsubjet = *it;
	}

	SDsubjetPt          = -std::numeric_limits<double>::max();
	SDsubjetEta         = -std::numeric_limits<double>::max();
	SDsubjetPhi         = -std::numeric_limits<double>::max();
	SDsubjetMass        = -std::numeric_limits<double>::max();
	SDsubjetDeepCSVb    = -std::numeric_limits<double>::max();
	SDsubjetDeepCSVc    = -std::numeric_limits<double>::max();
	SDsubjetDeepCSVudsg = -std::numeric_limits<double>::max();
	SDsubjetDeepCSVbb   = -std::numeric_limits<double>::max();
	SDdeltaRsubjetJet   = std::numeric_limits<double>::max();
       
	SDsubjetPt           = corrsubjet.pt();
	SDsubjetEta          = corrsubjet.eta();
	SDsubjetPhi          = corrsubjet.phi();
	SDsubjetMass         = corrsubjet.mass();
	SDsubjetDeepCSVb     = corrsubjet.bDiscriminator(bDiscriminant); 
	SDsubjetDeepCSVc     = corrsubjet.bDiscriminator("pfDeepCSVJetTags:probc");
	SDsubjetDeepCSVudsg  = corrsubjet.bDiscriminator("pfDeepCSVJetTags:probudsg");
	SDsubjetDeepCSVbb    = corrsubjet.bDiscriminator("pfDeepCSVJetTags:probbb");
	SDsubjetHFlav        = corrsubjet.hadronFlavour();
	SDsubjetBTag         = selector->isJetTagged(corrsubjet, event, true, 0, true);
	SDdeltaRsubjetJet    = deltaR(corrak8.eta(), corrak8.phi(), SDsubjetEta, SDsubjetPhi);

	if(SDsubjetDeepCSVb + SDsubjetDeepCSVbb > 0.1522) nSDSubsDeepCSVL++;
	if(SDsubjetBTag > 0) nSDSubsDeepCSVMSF++;
	if(selector->isJetTagged(corrsubjet, event, true, 1, true) > 0) nSDSubsDeepCSVM_bSFup++;
	if(selector->isJetTagged(corrsubjet, event, true, 2, true) > 0) nSDSubsDeepCSVM_bSFdn++;
	if(selector->isJetTagged(corrsubjet, event, true, 3, true) > 0) nSDSubsDeepCSVM_lSFup++;
	if(selector->isJetTagged(corrsubjet, event, true, 4, true) > 0) nSDSubsDeepCSVM_lSFdn++;

	theJetAK8SDSubjetPt.push_back(SDsubjetPt);
	theJetAK8SDSubjetEta.push_back(SDsubjetEta);
	theJetAK8SDSubjetPhi.push_back(SDsubjetPhi);
	theJetAK8SDSubjetMass.push_back(SDsubjetMass);
	theJetAK8SDSubjetCSVb.push_back(SDsubjetDeepCSVb);
	theJetAK8SDSubjetCSVc.push_back(SDsubjetDeepCSVc);
	theJetAK8SDSubjetCSVudsg.push_back(SDsubjetDeepCSVudsg);
	theJetAK8SDSubjetCSVbb.push_back(SDsubjetDeepCSVbb);
	theJetAK8SDSubjetHFlav.push_back(SDsubjetHFlav);
	theJetAK8SDSubjetBTag.push_back(SDsubjetBTag);
	theJetAK8SDSubjetDR.push_back(SDdeltaRsubjetJet);
      }

      theJetAK8SDSubjetIndex.push_back(SDSubJetIndex);
      theJetAK8SDSubjetSize.push_back(nSDSubJets);

      // Get CHS Soft drop subjets for subjet b-tagging
      theJetAK8SDSubjetNDeepCSVL.push_back(nSDSubsDeepCSVL);
      theJetAK8SDSubjetNDeepCSVMSF.push_back(nSDSubsDeepCSVMSF);
      theJetAK8SDSubjetNDeepCSVM_bSFup.push_back(nSDSubsDeepCSVM_bSFup);
      theJetAK8SDSubjetNDeepCSVM_bSFdn.push_back(nSDSubsDeepCSVM_bSFdn);
      theJetAK8SDSubjetNDeepCSVM_lSFup.push_back(nSDSubsDeepCSVM_lSFup);
      theJetAK8SDSubjetNDeepCSVM_lSFdn.push_back(nSDSubsDeepCSVM_lSFdn);

      double puppicorr = 1.0;
      float genCorr  = 1.;
      float recoCorr = 1.;
      float totalWeight = 1.;
	
      genCorr =  puppisd_corrGEN->Eval(corrak8.pt());
      if(fabs(corrak8.eta()) <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval(corrak8.pt());
      else recoCorr = puppisd_corrRECO_for->Eval(corrak8.pt());
	
      puppicorr = genCorr * recoCorr;
    
      theSoftDrop = puppi_softdrop.M();
      double theSoftDropCorrected = theSoftDrop*puppicorr;

      double ptscale_sd = 1.0;
      double ptscale_sd_up = 1.0;
      double ptscale_sd_dn = 1.0;
      double unc_sd_up = 1.0;
      double unc_sd_dn = 1.0;

      if(isMc){
	double res = 10.1/theSoftDropCorrected;
	double factor_sd = 1;
	double uncert_sd = 0.20;
	double factor_sd_up = factor_sd + uncert_sd;
	double factor_sd_dn = factor_sd - uncert_sd;
	
	if (factor_sd>1) {
	  JERrand.SetSeed(abs(static_cast<int>(corrak8.phi()*1e4)));	   
	  ptscale_sd = 1 + JERrand.Gaus(0,res)*sqrt(factor_sd*factor_sd - 1.0);
	}
	if (factor_sd_up>1) {
	  JERrand.SetSeed(abs(static_cast<int>(corrak8.phi()*1e4)));	   
	  ptscale_sd_up = 1 + JERrand.Gaus(0,res)*sqrt(factor_sd_up*factor_sd_up - 1.0);
	}
	if (factor_sd_dn>1) {
	  JERrand.SetSeed(abs(static_cast<int>(corrak8.phi()*1e4)));	   
	  ptscale_sd_dn = 1 + JERrand.Gaus(0,res)*sqrt(factor_sd_dn*factor_sd_dn - 1.0);
	}
	unc_sd_up = 1.0094; // + sqrt(unc*unc + 0.023*0.023);
	unc_sd_dn = 0.9906; //1 - sqrt(unc*unc + 0.023*0.023);
      }
      
      int MaxProb = 10;
      double doubleB = corrak8.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");

      if (theSoftDropCorrected < 105 && theSoftDropCorrected > 65 && theNjettinessTau2/theNjettinessTau1 < 0.55) MaxProb = 3;
      else if (theSoftDropCorrected > 135 && theSoftDropCorrected < 135 && doubleB > 0.6) MaxProb = 2;
      else if (theSoftDropCorrected > 105 && theSoftDropCorrected < 210 && theNjettinessTau3/theNjettinessTau2 < 0.65) MaxProb = 1;
      else if (nSDSubsDeepCSVMSF > 0) MaxProb = 4;
      else if (nSDSubsDeepCSVMSF == 0) MaxProb = 0;
      else MaxProb = 10;

      maxProb.push_back(MaxProb);

      theJetAK8SJIndex.push_back(SDSubJetIndex);
      theJetAK8SJSize.push_back(nSubJets);
      theJetAK8SoftDropRaw.push_back(theSoftDrop);
      theJetAK8SoftDropCorr.push_back(theSoftDropCorrected);
      theJetAK8SoftDrop.push_back(theSoftDropCorrected*ptscale_sd);
      theJetAK8SoftDrop_JMSup.push_back(theSoftDropCorrected*ptscale_sd*unc_sd_up);
      theJetAK8SoftDrop_JMSdn.push_back(theSoftDropCorrected*ptscale_sd*unc_sd_dn);
      theJetAK8SoftDrop_JMRup.push_back(theSoftDropCorrected*ptscale_sd_up);
      theJetAK8SoftDrop_JMRdn.push_back(theSoftDropCorrected*ptscale_sd_dn);


    }

    SetValue("maxProb", maxProb);
    SetValue("theJetAK8Pt",     theJetAK8Pt);
    SetValue("theJetAK8Eta",    theJetAK8Eta);
    SetValue("theJetAK8Phi",    theJetAK8Phi);
    SetValue("theJetAK8Energy", theJetAK8Energy);
    SetValue("theJetAK8CSV",    theJetAK8CSV);
    SetValue("theJetAK8DoubleB",    theJetAK8DoubleB);
    SetValue("theJetAK8JetCharge", theJetAK8JetCharge);
    SetValue("theJetAK8GenPt",  theJetAK8GenPt);
    SetValue("theJetAK8GenDR",  theJetAK8GenDR);
    SetValue("theJetAK8GenMass",  theJetAK8GenMass);

    SetValue("theJetAK8CHSPt",     theJetAK8CHSPt);
    SetValue("theJetAK8CHSEta",    theJetAK8CHSEta);
    SetValue("theJetAK8CHSPhi",    theJetAK8CHSPhi);
    SetValue("theJetAK8CHSMass", theJetAK8CHSMass);
    SetValue("theJetAK8SoftDropRaw", theJetAK8SoftDropRaw);
    SetValue("theJetAK8SoftDropCorr", theJetAK8SoftDropCorr);
    SetValue("theJetAK8SoftDrop", theJetAK8SoftDrop);
    SetValue("theJetAK8SoftDrop_JMSup", theJetAK8SoftDrop_JMSup);
    SetValue("theJetAK8SoftDrop_JMSdn", theJetAK8SoftDrop_JMSdn);
    SetValue("theJetAK8SoftDrop_JMRup", theJetAK8SoftDrop_JMRup);
    SetValue("theJetAK8SoftDrop_JMRdn", theJetAK8SoftDrop_JMRdn);

    SetValue("theJetAK8CHSPrunedMass",   theJetAK8CHSPrunedMass);
    SetValue("theJetAK8CHSSoftDropMass", theJetAK8CHSSoftDropMass);
    
    SetValue("theJetAK8NjettinessTau1", theJetAK8NjettinessTau1);
    SetValue("theJetAK8NjettinessTau2", theJetAK8NjettinessTau2);
    SetValue("theJetAK8NjettinessTau3", theJetAK8NjettinessTau3);
    SetValue("theJetAK8CHSTau1", theJetAK8CHSTau1);
    SetValue("theJetAK8CHSTau2", theJetAK8CHSTau2);
    SetValue("theJetAK8CHSTau3", theJetAK8CHSTau3);

    SetValue("theJetAK8SoftDropn2b1",theJetAK8SoftDropn2b1);
    SetValue("theJetAK8SoftDropn2b2",theJetAK8SoftDropn2b2);
    SetValue("theJetAK8SoftDropn3b1",theJetAK8SoftDropn3b1);
    SetValue("theJetAK8SoftDropn3b2",theJetAK8SoftDropn3b2);

    SetValue("theJetAK8Mass",   theJetAK8Mass);
    SetValue("theJetAK8nDaughters", theJetAK8nDaughters);

    SetValue("theJetAK8SDSubjetPt",   theJetAK8SDSubjetPt);   
    SetValue("theJetAK8SDSubjetEta",  theJetAK8SDSubjetEta);  
    SetValue("theJetAK8SDSubjetPhi",  theJetAK8SDSubjetPhi);  
    SetValue("theJetAK8SDSubjetMass", theJetAK8SDSubjetMass); 
    SetValue("theJetAK8SDSubjetCSVb",  theJetAK8SDSubjetCSVb); 
    SetValue("theJetAK8SDSubjetCSVc",  theJetAK8SDSubjetCSVc);
    SetValue("theJetAK8SDSubjetCSVudsg",  theJetAK8SDSubjetCSVudsg);
    SetValue("theJetAK8SDSubjetCSVbb",  theJetAK8SDSubjetCSVbb);
    SetValue("theJetAK8SDSubjetHFlav", theJetAK8SDSubjetHFlav);
    SetValue("theJetAK8SDSubjetBTag",  theJetAK8SDSubjetBTag);  
    SetValue("theJetAK8SDSubjetDR",   theJetAK8SDSubjetDR);   
    SetValue("theJetAK8SDSubjetIndex",theJetAK8SDSubjetIndex);
    SetValue("theJetAK8SDSubjetSize", theJetAK8SDSubjetSize); 

    SetValue("theJetAK8SDSubjetNDeepCSVL",theJetAK8SDSubjetNDeepCSVL);
    SetValue("theJetAK8SDSubjetNDeepCSVMSF",theJetAK8SDSubjetNDeepCSVMSF);
    SetValue("theJetAK8SDSubjetNDeepCSVM_bSFup",theJetAK8SDSubjetNDeepCSVM_bSFup);
    SetValue("theJetAK8SDSubjetNDeepCSVM_bSFdn",theJetAK8SDSubjetNDeepCSVM_bSFdn);
    SetValue("theJetAK8SDSubjetNDeepCSVM_lSFup",theJetAK8SDSubjetNDeepCSVM_lSFup);
    SetValue("theJetAK8SDSubjetNDeepCSVM_lSFdn",theJetAK8SDSubjetNDeepCSVM_lSFdn);

    SetValue("theJetAK8SJIndex",theJetAK8SJIndex);
    SetValue("theJetAK8SJSize", theJetAK8SJSize); 


    //////////////// TRUE HADRONIC W/Z/H/Top decays //////////////////
    std::vector<int>    HadronicVHtID;
    std::vector<int>    HadronicVHtStatus;
    std::vector<double> HadronicVHtPt;
    std::vector<double> HadronicVHtEta;
    std::vector<double> HadronicVHtPhi;
    std::vector<double> HadronicVHtEnergy;
    std::vector<double> HadronicVHtD0Pt;
    std::vector<double> HadronicVHtD0Eta;
    std::vector<double> HadronicVHtD0Phi;
    std::vector<double> HadronicVHtD0E;
    std::vector<double> HadronicVHtD1Pt;
    std::vector<double> HadronicVHtD1Eta;
    std::vector<double> HadronicVHtD1Phi;
    std::vector<double> HadronicVHtD1E;
    std::vector<double> HadronicVHtD2Pt;
    std::vector<double> HadronicVHtD2Eta;
    std::vector<double> HadronicVHtD2Phi;
    std::vector<double> HadronicVHtD2E;
  
    // Get the generated particle collection
    edm::Handle<reco::GenParticleCollection> genParticles;
    if(isMc && event.getByLabel(genParticles_it, genParticles)){
      
      for(size_t i = 0; i < genParticles->size(); i++){
	const reco::GenParticle &p = (*genParticles).at(i);
	int id = p.pdgId();

	bool hasRadiation = false;
	bool hasLepton = false;
	  
	if(abs(id) == 23 || abs(id) == 24 || abs(id) == 25 || abs(id) == 6){
	    
	  size_t nDs = p.numberOfDaughters();
	  for(size_t j = 0; j < nDs; j++){
	    int dauId = (p.daughter(j))->pdgId();
	    const reco::Candidate *d = p.daughter(j);
	    if(d->pdgId() != dauId) std::cout << "making daughter GenParticle didn't work" << std::endl;
	    
	    if(abs(dauId) == abs(id)) hasRadiation = true;
	    else if(abs(dauId) == 24){  // check t->Wb->leptons and H->WW->leptons
	      while(d->numberOfDaughters() == 1) d = d->daughter(0);
	      if(abs(d->daughter(0)->pdgId()) > 10 && abs(d->daughter(0)->pdgId()) < 17) hasLepton = true;
	      if(abs(d->daughter(1)->pdgId()) > 10 && abs(d->daughter(1)->pdgId()) < 17) hasLepton = true;
	    }
	    else if(abs(dauId) == 23){  // check H->ZZ->leptons
	      while(d->numberOfDaughters() == 1) d = d->daughter(0);
	      if(abs(d->daughter(0)->pdgId()) > 10 && abs(d->daughter(0)->pdgId()) < 17) hasLepton = true;
	      if(abs(d->daughter(1)->pdgId()) > 10 && abs(d->daughter(1)->pdgId()) < 17) hasLepton = true;
	    }
	    else if(abs(dauId) > 10 && abs(dauId) < 17) hasLepton = true;	    
	    
	  }
	  
	  if(hasRadiation) continue;	  
	  if(hasLepton) continue;	  
	  if(p.pt() < 175) continue;
	  
	  if(abs(id) == 24){
	    double dRWb = 1000;
	    double dRWW = 1000;
	    
	    const reco::Candidate *mother = p.mother();
	    while(abs(mother->pdgId()) == 24) mother = mother->mother();
	    
	    if(abs(mother->pdgId()) == 6){
	      double dr = reco::deltaR(p.p4(),mother->daughter(1)->p4());
	      if(abs(mother->daughter(1)->pdgId()) == 24) dr = reco::deltaR(p.p4(),mother->daughter(0)->p4());
	      if(dr < dRWb) dRWb = dr;
	    }else if(abs(mother->pdgId()) == 25){
	      double dr = 1000;
	      if(p.pdgId()*mother->daughter(0)->pdgId() > 0){
		dr = reco::deltaR(p.p4(),mother->daughter(1)->p4());
	      }else{
		dr = reco::deltaR(p.p4(),mother->daughter(0)->p4());
	      }
	      if(dr < dRWW) dRWW = dr;
	    }		
	    
	    if(dRWW < 0.8) continue; // W from merged H
	    if(dRWb < 0.8) continue; // W from merged t
	  }
	  
	  if(abs(id) == 23){
	    double dRZZ = 1000;
	    
	    const reco::Candidate *mother = p.mother();
	    while(abs(mother->pdgId()) == 23) mother = mother->mother();
	    
	    if(abs(mother->pdgId()) == 25){
	      double dr = 1000;
	      if(p.pdgId()*mother->daughter(0)->pdgId() > 0){
		dr = reco::deltaR(p.p4(),mother->daughter(1)->p4());
	      }else{
		dr = reco::deltaR(p.p4(),mother->daughter(0)->p4());
	      }
	      if(dr < dRZZ) dRZZ = dr;
	    }
	    
	    if(dRZZ < 0.8) continue; // Z from merged H
	  }
	  
	  if(p.numberOfDaughters() < 2){
	    std::cout << p.numberOfDaughters() << " daughters from " << p.pdgId() << std::endl;
	    continue;
	  }

	  HadronicVHtStatus.push_back( p.status() );
	  HadronicVHtID.push_back( p.pdgId() );
	  HadronicVHtPt.push_back( p.pt() );
	  HadronicVHtEta.push_back( p.eta() );
	  HadronicVHtPhi.push_back( p.phi() );
	  HadronicVHtEnergy.push_back( p.energy() );
	  
	  if(abs(id) != 6){
	    HadronicVHtD0Pt.push_back( p.daughter(0)->pt());
	    HadronicVHtD0Eta.push_back( p.daughter(0)->eta());
	    HadronicVHtD0Phi.push_back( p.daughter(0)->phi());
	    HadronicVHtD0E.push_back( p.daughter(0)->energy());
	    HadronicVHtD1Pt.push_back( p.daughter(1)->pt());	    
	    HadronicVHtD1Eta.push_back( p.daughter(1)->eta());
	    HadronicVHtD1Phi.push_back( p.daughter(1)->phi());
	    HadronicVHtD1E.push_back( p.daughter(1)->energy());
	    HadronicVHtD2Pt.push_back(-99.9);
	    HadronicVHtD2Eta.push_back(-99.9);
	    HadronicVHtD2Phi.push_back(-99.9);
	    HadronicVHtD2E.push_back(-99.9);
	  }else{
	    const reco::Candidate *W = p.daughter(0);
	    const reco::Candidate *b = p.daughter(1);
	    if(fabs(W->pdgId()) != 24){
	      W = p.daughter(1);
	      b = p.daughter(0);
	    }
	    while(W->numberOfDaughters() == 1) W = W->daughter(0);
	    if(W->daughter(1)->pdgId() == 22) W = W->daughter(0);	    
	    while(W->numberOfDaughters() == 1) W = W->daughter(0);
	    if(W->daughter(1)->pdgId() == 22) W = W->daughter(0);
            while(W->numberOfDaughters() == 1) W = W->daughter(0);
	    if(W->daughter(1)->pdgId()==22) cout << "weird W decay to photons" << endl;

	    HadronicVHtD0Pt.push_back( b->pt());
	    HadronicVHtD0Eta.push_back( b->eta());
	    HadronicVHtD0Phi.push_back( b->phi());
	    HadronicVHtD0E.push_back( b->energy());
	    HadronicVHtD1Pt.push_back( W->daughter(0)->pt());
	    HadronicVHtD2Pt.push_back( W->daughter(1)->pt());
	    HadronicVHtD1Eta.push_back( W->daughter(0)->eta());
	    HadronicVHtD2Eta.push_back( W->daughter(1)->eta());
	    HadronicVHtD1Phi.push_back( W->daughter(0)->phi());
	    HadronicVHtD2Phi.push_back( W->daughter(1)->phi());
	    HadronicVHtD1E.push_back( W->daughter(0)->energy());
	    HadronicVHtD2E.push_back( W->daughter(1)->energy());
	  }
	}
      }
    }

    SetValue("HadronicVHtStatus",HadronicVHtStatus);
    SetValue("HadronicVHtID",HadronicVHtID);
    SetValue("HadronicVHtPt",HadronicVHtPt);
    SetValue("HadronicVHtEta",HadronicVHtEta);
    SetValue("HadronicVHtPhi",HadronicVHtPhi);
    SetValue("HadronicVHtEnergy",HadronicVHtEnergy);
    SetValue("HadronicVHtD0Pt",HadronicVHtD0Pt);
    SetValue("HadronicVHtD0Eta",HadronicVHtD0Eta);
    SetValue("HadronicVHtD0Phi",HadronicVHtD0Phi);
    SetValue("HadronicVHtD0E",HadronicVHtD0E);
    SetValue("HadronicVHtD1Pt",HadronicVHtD1Pt);
    SetValue("HadronicVHtD1Eta",HadronicVHtD1Eta);
    SetValue("HadronicVHtD1Phi",HadronicVHtD1Phi);
    SetValue("HadronicVHtD1E",HadronicVHtD1E);
    SetValue("HadronicVHtD2Pt",HadronicVHtD2Pt);
    SetValue("HadronicVHtD2Eta",HadronicVHtD2Eta);
    SetValue("HadronicVHtD2Phi",HadronicVHtD2Phi);
    SetValue("HadronicVHtD2E",HadronicVHtD2E);

    return 0;
}

int JetSubCalc::EndJob()
{
  delete jecUnc;
  delete jecak8;
  delete L3JetParAK8MC;
  delete L2JetParAK8MC;
  delete ResJetParAK8MC; 
  delete L3JetParAK8BCD;
  delete L2JetParAK8BCD;
  delete ResJetParAK8BCD; 
  delete L3JetParAK8EF;
  delete L2JetParAK8EF;
  delete ResJetParAK8EF; 
  delete L3JetParAK8G;
  delete L2JetParAK8G;
  delete ResJetParAK8G; 
  delete L3JetParAK8H;
  delete L2JetParAK8H;
  delete ResJetParAK8H; 
  return 0;
}
