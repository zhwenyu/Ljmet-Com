/*
 Calculator for substructure variables
 
 Author: Joshua Swanson, 2014
 */

#include <iostream>
#include <limits>   // std::numeric_limits

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"

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
    std::string bDiscriminant;
    std::string tagInfo;
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
    if (mPset.exists("slimmedJetColl")) slimmedJetColl_it = mPset.getParameter<edm::InputTag>("slimmedJetColl");
    else slimmedJetColl_it = edm::InputTag("slimmedJet");
    
    if (mPset.exists("slimmedJetsAK8Coll")) slimmedJetsAK8Coll_it = mPset.getParameter<edm::InputTag>("slimmedJetsAK8Coll");
    else slimmedJetsAK8Coll_it = edm::InputTag("slimmedJetsAK8");
    
    if (mPset.exists("bDiscriminant")) bDiscriminant = mPset.getParameter<std::string>("bDiscriminant");
    else bDiscriminant = "combinedSecondaryVertexBJetTags";
    
    if (mPset.exists("tagInfo")) tagInfo = mPset.getParameter<std::string>("tagInfo");
    else tagInfo = "CATop";
    
    return 0;
}

int JetSubCalc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector)
{
    edm::Handle<std::vector<pat::Jet> > theAK8Jets;
    event.getByLabel(slimmedJetsAK8Coll_it, theAK8Jets);
    
    // Four vector
    std::vector<double> theJetAK8Pt;
    std::vector<double> theJetAK8Eta;
    std::vector<double> theJetAK8Phi;
    std::vector<double> theJetAK8Energy;
    
    std::vector<double> theJetAK8Mass;
    
    std::vector<int>    theJetAK8Index;
    std::vector<double> theJetAK8CSV;
    std::vector<int>    theJetAK8nDaughters;
    
    // Daughter four vector and index
    std::vector<double> theJetAK8DaughterPt;
    std::vector<double> theJetAK8DaughterEta;
    std::vector<double> theJetAK8DaughterPhi;
    std::vector<double> theJetAK8DaughterEnergy;
    
    std::vector<int> theJetAK8DaughterMotherIndex;
    
    std::vector<int> theJetAK8CSVLSubJets;
    std::vector<int> theJetAK8CSVMSubJets;
    std::vector<int> theJetAK8CSVTSubJets;
    
    float subjetCSV;
    int CSVL, CSVM, CSVT;
    
    for (std::vector<pat::Jet>::const_iterator ijet = theAK8Jets->begin(); ijet != theAK8Jets->end(); ijet++) {
        int index = (int)(ijet-theAK8Jets->begin());
        
        subjetCSV = -std::numeric_limits<float>::max();
        CSVL = 0;
        CSVM = 0;
        CSVT = 0;
        
        theJetAK8Pt.push_back(ijet->pt());
        theJetAK8Eta.push_back(ijet->eta());
        theJetAK8Phi.push_back(ijet->phi());
        theJetAK8Energy.push_back(ijet->energy());
        
        theJetAK8Mass.push_back(ijet->mass());
        theJetAK8CSV.push_back(ijet->bDiscriminator( bDiscriminant ));
        theJetAK8nDaughters.push_back((int)ijet->numberOfDaughters());
        
        theJetAK8Index.push_back(index);
        
        for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++) {
            theJetAK8DaughterPt.push_back(ijet->daughter(ui)->pt());
            theJetAK8DaughterEta.push_back(ijet->daughter(ui)->eta());
            theJetAK8DaughterPhi.push_back(ijet->daughter(ui)->phi());
            theJetAK8DaughterEnergy.push_back(ijet->daughter(ui)->energy());
            
            theJetAK8DaughterMotherIndex.push_back(index);
            
            pat::PackedCandidate const * subjet = dynamic_cast<pat::PackedCandidate const *>(ijet->daughter(ui));
//            subjetCSV = subjet->bDiscriminator(bDiscriminant);
//            if (subjetCSV > 0.244 && ijet->daughter(ui)->pt() > 20.) {
//                CSVL++;
//            }
//            if (subjetCSV > 0.679 && ijet->daughter(ui)->pt() > 20.){
//                CSVM++;
//            }
//            if (subjetCSV > 0.898 && ijet->daughter(ui)->pt() > 20.){
//                CSVT++;
//            }
        }
        theJetAK8CSVLSubJets.push_back(CSVL);
        theJetAK8CSVMSubJets.push_back(CSVM);
        theJetAK8CSVTSubJets.push_back(CSVT);
    }
    
    // Four vector
    SetValue("theJetAK8Pt",     theJetAK8Pt);
    SetValue("theJetAK8Eta",    theJetAK8Eta);
    SetValue("theJetAK8Phi",    theJetAK8Phi);
    SetValue("theJetAK8Energy", theJetAK8Energy);
    
    SetValue("theJetAK8Mass",   theJetAK8Mass);
    
    SetValue("theJetAK8Index",      theJetAK8Index);
    SetValue("theJetAK8CSV",        theJetAK8CSV);
    SetValue("theJetAK8nDaughters", theJetAK8nDaughters);
    
    // Daughter four vector and index
    SetValue("theJetAK8DaughterPt",     theJetAK8DaughterPt);
    SetValue("theJetAK8DaughterEta",    theJetAK8DaughterEta);
    SetValue("theJetAK8DaughterPhi",    theJetAK8DaughterPhi);
    SetValue("theJetAK8DaughterEnergy", theJetAK8DaughterEnergy);
    
    SetValue("theJetAK8DaughterMotherIndex", theJetAK8DaughterMotherIndex);
    
    SetValue("theJetAK8CSVLSubJets", theJetAK8CSVLSubJets);
    SetValue("theJetAK8CSVMSubJets", theJetAK8CSVMSubJets);
    SetValue("theJetAK8CSVTSubJets", theJetAK8CSVTSubJets);
    
    edm::Handle<std::vector<pat::Jet> > theJets;
    event.getByLabel(slimmedJetColl_it, theJets);
    
    //Four vector
    std::vector<double> theJetPt;
    std::vector<double> theJetEta;
    std::vector<double> theJetPhi;
    std::vector<double> theJetEnergy;
    
    std::vector<double> theJetCSV;
    //std::vector<double> theJetRCN;
    
    //Identity
    std::vector<int> theJetIndex;
    std::vector<int> theJetnDaughters;
    
    //Top-like properties
    std::vector<double> theJetTopMass;
    std::vector<double> theJetMinPairMass;
    std::vector<bool> theJetIsTopTagged;
    
    //Daughter four vector and index
    std::vector<double> theJetDaughterPt;
    std::vector<double> theJetDaughterEta;
    std::vector<double> theJetDaughterPhi;
    std::vector<double> theJetDaughterEnergy;
    
    std::vector<int> theJetDaughterMotherIndex;
    
    std::vector<int> theJetCSVLSubJets;
    std::vector<int> theJetCSVMSubJets;
    std::vector<int> theJetCSVTSubJets;

    bool topTagged;
    double topMass, minMass;
    int nSubJets;

    for (std::vector<pat::Jet>::const_iterator ijet = theJets->begin(); ijet != theJets->end(); ijet++){
        int index = (int)(ijet-theJets->begin());
        
        subjetCSV = -std::numeric_limits<float>::max();
        CSVL = 0;
        CSVM = 0;
        CSVT = 0;
        
        theJetPt.push_back(ijet->pt());
        theJetEta.push_back(ijet->eta());
        theJetPhi.push_back(ijet->phi());
        theJetEnergy.push_back(ijet->energy());
        
        theJetCSV.push_back(ijet->bDiscriminator( bDiscriminant ));
        //theJetRCN.push_back((ijet->chargedEmEnergy()+ijet->chargedHadronEnergy()) / (ijet->neutralEmEnergy()+ijet->neutralHadronEnergy()));
        
        theJetIndex.push_back(index);
        theJetnDaughters.push_back((int)ijet->numberOfDaughters());
        
//        reco::CATopJetTagInfo* jetInfo = (reco::CATopJetTagInfo*) ijet->tagInfo( tagInfo );
        reco::CATopJetTagInfo const * jetInfo = dynamic_cast<reco::CATopJetTagInfo const *>( ijet->tagInfo("caTop"));
        
        topTagged = false;
        topMass   = -std::numeric_limits<double>::max();
        minMass   = -std::numeric_limits<double>::max();
        nSubJets  = std::numeric_limits<int>::min();
        
        if ( jetInfo != 0 ) {
            topMass = jetInfo->properties().topMass;
            minMass = jetInfo->properties().minMass;
            nSubJets = jetInfo->properties().nSubJets;
            std::cout << topMass << '\t' << minMass << '\t' << nSubJets << std::endl;
            
            if ( nSubJets > 2 && minMass > 50.0 && topMass > 150.0 ) {
                topTagged = true;
            }
        }
        
        theJetTopMass.push_back(topMass);
        theJetMinPairMass.push_back(minMass);
        theJetIsTopTagged.push_back(topTagged);
        
        for (size_t ui = 0; ui < ijet->numberOfDaughters(); ui++){
            theJetDaughterPt . push_back(ijet->daughter(ui)->pt());
            theJetDaughterEta . push_back(ijet->daughter(ui)->eta());
            theJetDaughterPhi . push_back(ijet->daughter(ui)->phi());
            theJetDaughterEnergy . push_back(ijet->daughter(ui)->energy());
            
            theJetDaughterMotherIndex . push_back(index);
            pat::PackedCandidate const * subjet = dynamic_cast<pat::PackedCandidate const *>(ijet->daughter(ui));

//            pat::Jet const * subjet = dynamic_cast<pat::Jet const *>(ijet->daughter(ui));
//            subjetCSV = subjet->bDiscriminator(bDiscriminant);
//            if (subjetCSV > 0.244 && ijet->daughter(ui)->pt() > 20.) {
//                CSVL++;
//            }
//            if (subjetCSV > 0.679 && ijet->daughter(ui)->pt() > 20.) {
//                CSVM++;
//            }
//            if (subjetCSV > 0.898 && ijet->daughter(ui)->pt() > 20.) {
//                CSVT++;
//            }
        }
        theJetCSVLSubJets.push_back(CSVL);
        theJetCSVMSubJets.push_back(CSVM);
        theJetCSVTSubJets.push_back(CSVT);
    }
    
    // Four vector
    SetValue("theJetPt",     theJetPt);
    SetValue("theJetEta",    theJetEta);
    SetValue("theJetPhi",    theJetPhi);
    SetValue("theJetEnergy", theJetEnergy);
    
    SetValue("theJetCSV",    theJetCSV);
    //SetValue("theJetRCN",  theJetRCN);
    
    // Identity
    SetValue("theJetIndex",      theJetIndex);
    SetValue("theJetnDaughters", theJetnDaughters);
    
    // Properties
    SetValue("theJetTopMass",     theJetTopMass);
    SetValue("theJetMinPairMass", theJetMinPairMass);
    SetValue("theJetIsTopTagged", theJetIsTopTagged);
    
    // Daughter four vector and index
    SetValue("theJetDaughterPt",     theJetDaughterPt);
    SetValue("theJetDaughterEta",    theJetDaughterEta);
    SetValue("theJetDaughterPhi",    theJetDaughterPhi);
    SetValue("theJetDaughterEnergy", theJetDaughterEnergy);
    
    SetValue("theJetDaughterMotherIndex", theJetDaughterMotherIndex);
    
    SetValue("theJetCSVLSubJets", theJetCSVLSubJets);
    SetValue("theJetCSVMSubJets", theJetCSVMSubJets);
    SetValue("theJetCSVTSubJets", theJetCSVTSubJets);
    
    return 0;
}

int JetSubCalc::EndJob()
{
    return 0;
}