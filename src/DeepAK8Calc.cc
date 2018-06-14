#include <iostream>
#include <algorithm>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "TLorentzVector.h" 
#include "DataFormats/Math/interface/deltaR.h"
#include "LJMet/Com/interface/DeepAK8Calc.h"

static int reg = LjmetFactory::GetInstance()->Register(new DeepAK8Calc(), "DeepAK8Calc");

using namespace std;

DeepAK8Calc::DeepAK8Calc(){

}

int DeepAK8Calc::BeginJob(){

  return 0;

}


int DeepAK8Calc::AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector){
  //DECLARE VARIABLES
  std::vector<double> dnn_B;
  std::vector<double> dnn_J;
  std::vector<double> dnn_W;
  std::vector<double> dnn_Z;
  std::vector<double> dnn_H;
  std::vector<double> dnn_T;
  
  std::vector<double> decorr_B;
  std::vector<double> decorr_J;
  std::vector<double> decorr_W;
  std::vector<double> decorr_Z;
  std::vector<double> decorr_H;
  std::vector<double> decorr_T;

  //getByLabel
  edm::InputTag AK8JetColl = edm::InputTag("slimmedJetsAK8");
  edm::Handle<std::vector<pat::Jet> > AK8Jets;
  event.getByLabel(AK8JetColl, AK8Jets);

  edm::InputTag DeepAK8rawB = edm::InputTag("deepntuplizer:DeepAK8rawB");
  edm::Handle<std::vector<float> > DeepAK8_raw_B;
  event.getByLabel(DeepAK8rawB, DeepAK8_raw_B);

  edm::InputTag DeepAK8rawJ= edm::InputTag("deepntuplizer:DeepAK8rawJ");
  edm::Handle<std::vector<float> > DeepAK8_raw_J;
  event.getByLabel(DeepAK8rawJ, DeepAK8_raw_J);

  edm::InputTag DeepAK8rawW = edm::InputTag("deepntuplizer:DeepAK8rawW");
  edm::Handle<std::vector<float> > DeepAK8_raw_W;
  event.getByLabel(DeepAK8rawW, DeepAK8_raw_W);

  edm::InputTag DeepAK8rawZ = edm::InputTag("deepntuplizer:DeepAK8rawZ");
  edm::Handle<std::vector<float> > DeepAK8_raw_Z;
  event.getByLabel(DeepAK8rawZ, DeepAK8_raw_Z);

  edm::InputTag DeepAK8rawH = edm::InputTag("deepntuplizer:DeepAK8rawH");
  edm::Handle<std::vector<float> > DeepAK8_raw_H;
  event.getByLabel(DeepAK8rawH, DeepAK8_raw_H);

  edm::InputTag DeepAK8rawT = edm::InputTag("deepntuplizer:DeepAK8rawT");
  edm::Handle<std::vector<float> > DeepAK8_raw_T;
  event.getByLabel(DeepAK8rawT, DeepAK8_raw_T);

  edm::InputTag DeepAK8DecorrrawB = edm::InputTag("deepntuplizer:DeepAK8DecorrrawB");
  edm::Handle<std::vector<float> > DeepAK8Decorr_raw_B;
  event.getByLabel(DeepAK8DecorrrawB, DeepAK8Decorr_raw_B);

  edm::InputTag DeepAK8DecorrrawJ = edm::InputTag("deepntuplizer:DeepAK8DecorrrawJ");
  edm::Handle<std::vector<float> > DeepAK8Decorr_raw_J;
  event.getByLabel(DeepAK8DecorrrawJ, DeepAK8Decorr_raw_J);

  edm::InputTag DeepAK8DecorrrawW = edm::InputTag("deepntuplizer:DeepAK8DecorrrawW");
  edm::Handle<std::vector<float> > DeepAK8Decorr_raw_W;
  event.getByLabel(DeepAK8DecorrrawW, DeepAK8Decorr_raw_W);

  edm::InputTag DeepAK8DecorrrawZ = edm::InputTag("deepntuplizer:DeepAK8DecorrrawZ");
  edm::Handle<std::vector<float> > DeepAK8Decorr_raw_Z;
  event.getByLabel(DeepAK8DecorrrawZ, DeepAK8Decorr_raw_Z);

  edm::InputTag DeepAK8DecorrrawH = edm::InputTag("deepntuplizer:DeepAK8DecorrrawH");
  edm::Handle<std::vector<float> > DeepAK8Decorr_raw_H;
  event.getByLabel(DeepAK8DecorrrawH, DeepAK8Decorr_raw_H);

  edm::InputTag DeepAK8DecorrrawT = edm::InputTag("deepntuplizer:DeepAK8DecorrrawT");
  edm::Handle<std::vector<float> > DeepAK8Decorr_raw_T;
  event.getByLabel(DeepAK8DecorrrawT, DeepAK8Decorr_raw_T);

  //Size Check to make sure AK8Jets and the b, j, W... vectors are the same length
  if(DeepAK8_raw_B->size() != AK8Jets->size()) {
    cerr << "WARNING: THE DeepAK8Jets VECTOR IS SIZE " << AK8Jets->size() << ", AND THE OTHER VECTORS ARE SIZE " << DeepAK8_raw_B->size() << endl;
  }

  int i = 0; //i is an iterator
  for (std::vector<pat::Jet>::const_iterator ijet = AK8Jets->begin(); ijet != AK8Jets->end(); ijet++){
    
    if(ijet->pt() < 200) { // not all info there for lower pt
      i++;
      continue;
    }

    //PF Tight    
    bool tightJetID = false;
    pat::Jet rawJet = ijet->correctedJet(0);
    if(abs(rawJet.eta()) <= 2.7){
      tightJetID = (rawJet.neutralHadronEnergyFraction() < 0.9 &&
                    rawJet.neutralEmEnergyFraction() < 0.9 &&
                    (rawJet.chargedMultiplicity()+rawJet.neutralMultiplicity()) > 1) &&
        ((abs(rawJet.eta()) <= 2.4 &&
          rawJet.chargedHadronEnergyFraction() > 0 &&
          //rawJet.chargedEmEnergyFraction() < 0.99 &&                                                                                                                                           
          rawJet.chargedMultiplicity() > 0) ||
         abs(rawJet.eta()) > 2.4);
    }else{
      tightJetID = true;
    }

    if(!tightJetID) {
      i++;
      continue;
    }

    //PUSH_BACK(...)
    dnn_B.push_back(DeepAK8_raw_B->at(i));
    dnn_J.push_back(DeepAK8_raw_J->at(i));
    dnn_W.push_back(DeepAK8_raw_W->at(i));
    dnn_Z.push_back(DeepAK8_raw_Z->at(i));
    dnn_H.push_back(DeepAK8_raw_H->at(i));
    dnn_T.push_back(DeepAK8_raw_T->at(i));

    decorr_B.push_back(DeepAK8Decorr_raw_B->at(i));
    decorr_J.push_back(DeepAK8Decorr_raw_J->at(i));
    decorr_W.push_back(DeepAK8Decorr_raw_W->at(i));
    decorr_Z.push_back(DeepAK8Decorr_raw_Z->at(i));
    decorr_H.push_back(DeepAK8Decorr_raw_H->at(i));
    decorr_T.push_back(DeepAK8Decorr_raw_T->at(i));

    i++;
  }

  //SETVALUES...
  SetValue("dnn_B",dnn_B);
  SetValue("dnn_J",dnn_J);
  SetValue("dnn_W",dnn_W);
  SetValue("dnn_Z",dnn_Z);
  SetValue("dnn_H",dnn_H);
  SetValue("dnn_T",dnn_T);

  SetValue("decorr_B",decorr_B);
  SetValue("decorr_J",decorr_J);
  SetValue("decorr_W",decorr_W);
  SetValue("decorr_Z",decorr_Z);
  SetValue("decorr_H",decorr_H);
  SetValue("decorr_T",decorr_T);

  return 0;

}


int DeepAK8Calc::EndJob()
{
  return 0;
}
