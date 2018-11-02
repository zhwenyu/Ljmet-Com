#ifndef LJMet_Com_interface_DeepAK8Calc_h
#define LJMet_Com_interface_DeepAK8Calc_h

/*
 Best class for all calculators
 
 Author: Gena Kukartsev, 2012
 */

#include <iostream>
#include <vector>

#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
//#include "LJMet/Com/interface/VVString.h"

// CMS
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"


using namespace std;

class LjmetFactory;

class DeepAK8Calc : public BaseCalc {
    //
    // DeepAK8 class for all calculators
    //
    //

    
 public:
    DeepAK8Calc();
    virtual ~DeepAK8Calc() { }
    DeepAK8Calc(const DeepAK8Calc &); // stop default
    std::string GetName() { return mName; }
    virtual int BeginJob();
    virtual int ProduceEvent(edm::EventBase const & event, BaseEventSelector * selector) { return 0; }
    virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
    virtual int EndJob();

 private:

};

#endif
