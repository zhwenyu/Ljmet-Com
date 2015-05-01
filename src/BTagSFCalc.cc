/*
 Calculator for many btagging algorithms
 
 Author: Alex Garabedian, 2014
 */

#include <iostream>
#include <map>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

class LjmetFactory;

class BTagSFCalc : public BaseCalc{
public:
    
    BTagSFCalc();
    virtual ~BTagSFCalc(){}
    
    virtual int BeginJob(){
        
        if (mPset.exists("triggerSummary")) triggerSummary_ = mPset.getParameter<edm::InputTag>("triggerSummary");
        else                                triggerSummary_ = edm::InputTag("hltTriggerSummaryAOD");
        
        if (mPset.exists("rhoSrc")) rhoSrc_ = mPset.getParameter<edm::InputTag>("rhoSrc");
        else                        rhoSrc_ = edm::InputTag("kt6PFJets", "rho");
        
        if (mPset.exists("isTB"))    isTB_ = mPset.getParameter<bool>("isTB");
        else                         isTB_ = false;
        
        if (mPset.exists("isTT"))    isTT_ = mPset.getParameter<bool>("isTT");
        else                         isTT_ = false;
        
        if (mPset.exists("isWJets")) isWJets_ = mPset.getParameter<bool>("isWJets");
        else                         isWJets_ = false;
        
        return 0;
    }
    
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob(){return 0;}
    
    
private:
    
    edm::InputTag triggerSummary_;
    edm::InputTag rhoSrc_;
    bool isWJets_;
    bool isTB_;
    bool isTT_;
    
};

static int reg = LjmetFactory::GetInstance()->Register(new BTagSFCalc(), "BTagSFCalc");

BTagSFCalc::BTagSFCalc(){
}

int BTagSFCalc::AnalyzeEvent(edm::EventBase const & event,
                                BaseEventSelector * selector){
    
    int _nJets = 8;
    //
    // compute event variables here
    //
    
    //
    // _____ Get objects from the selector _____________________
    std::vector<edm::Ptr<pat::Jet> >            const & vSelJets = selector->GetSelectedJets();
    std::vector<edm::Ptr<pat::Jet> >            const & vAllJets = selector->GetAllJets();
    
    string tagger_tche = "trackCountingHighEffBJetTags";
    string tagger_tchp = "trackCountingHighPurBJetTags";
    string tagger_ssv = "simpleSecondaryVertexBJetTags";
    string tagger_ssv_hp = "simpleSecondaryVertexHighPurBJetTags";
    string tagger_ssv_he = "simpleSecondaryVertexHighEffBJetTags";
    string tagger_jBp = "jetBProbabilityBJetTags";
    string tagger_jp = "jetProbabilityBJetTags";
    string tagger_csv = "combinedSecondaryVertexV2BJetTags";
    string tagger_csvmva = "combinedSecondaryVertexV2MVABJetTags";
    string tagger_csvivf = "combinedInclusiveSecondaryVertexV2BJetTags";
    string tagger_cmva = "combinedMVABJetTags";
    
    vector<string> tagger;
    tagger.push_back(tagger_tche);
    tagger.push_back(tagger_tchp);
    tagger.push_back(tagger_ssv);
    tagger.push_back(tagger_ssv_hp);
    tagger.push_back(tagger_ssv_he);
    tagger.push_back(tagger_jBp);
    tagger.push_back(tagger_jp);
    tagger.push_back(tagger_csv);
    tagger.push_back(tagger_csvmva);
    tagger.push_back(tagger_csvivf);
    tagger.push_back(tagger_cmva);
    
    
    map< string, vector<Double_t> > Discrims;
    
    for( Int_t jetcount = 0; jetcount < _nJets; jetcount++){
        for(vector<string>::iterator st  = tagger.begin(); st != tagger.end(); ++st){
            Discrims[*st].push_back(-15.0);
        }
    }
    
    int jetNum=0;
    for(vector<edm::Ptr<pat::Jet> >::const_iterator jet = vSelJets.begin(); jet != vSelJets.end(); ++jet){
        for(vector<string>::iterator st  = tagger.begin(); st != tagger.end(); ++st){
            Discrims[*st][jetNum] = (*jet)->bDiscriminator(*st);
        }
        if (jetNum == (_nJets - 1)) break;
        jetNum++;
    }
    
    for(vector<string>::iterator st  = tagger.begin(); st != tagger.end(); ++st){
        SetValue(*st, Discrims[*st]);
    }
    
    
    
    
    
    return 0;
}
