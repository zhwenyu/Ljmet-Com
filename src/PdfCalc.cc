/*
  Calculator for PDF weights and uncertainties

   Author: Gena Kukartsev, 2013
*/



#include <iostream>
#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetFactory.h"
#include "LJMet/Com/interface/LjmetEventContent.h"
#include "LJMet/Com/interface/LjmetPdfWeightProducer.h"



class LjmetFactory;



class PdfCalc : public BaseCalc{
  
 public:
  
  PdfCalc();
  virtual ~PdfCalc();

  virtual int BeginJob();
  virtual int AnalyzeEvent(edm::EventBase const & event, BaseEventSelector * selector);
  virtual int EndJob();

  
 private:
  
  LjmetPdfWeightProducer * pPdfWeights;


};



static int reg = LjmetFactory::GetInstance()->Register(new PdfCalc(), "PdfCalc");



PdfCalc::PdfCalc(){
}



PdfCalc::~PdfCalc(){
}



int PdfCalc::BeginJob(){

  // initialize
  pPdfWeights = new LjmetPdfWeightProducer(mPset);
  pPdfWeights->beginJob();

  return 0;
}



int PdfCalc::EndJob(){

  delete pPdfWeights;

  return 0;
}



int PdfCalc::AnalyzeEvent(edm::EventBase const & event,
			     BaseEventSelector * selector){
  //
  // compute event variables here
  //


  // use PDF weight producer to calculate vectors of weights
  std::map<std::string,std::vector<double> > mPdfs = pPdfWeights->produce(event);


  // for each PDF set, save weights, averages etc.
  std::map<std::string,std::vector<double> >::const_iterator iPdf;
  for (iPdf=mPdfs.begin(); iPdf!=mPdfs.end(); ++iPdf){
    std::string pdfName = iPdf->first;
    std::vector<double> const & vWeights = iPdf->second;

    // save the vector of weights
    SetValue("PdfWeightsVec_"+pdfName, vWeights);

    // compute and save average +- weights
    unsigned int nPdfs = vWeights.size();
    double weight_average = vWeights[0];
    double weight_plus = 0;
    double weight_minus = 0;
    if (nPdfs>=3){
      for (unsigned int i=1; i!=nPdfs; i+=2){
	weight_plus  += vWeights[i];
	weight_minus += vWeights[i+1];
	weight_average += vWeights[i] + vWeights[i+1];
      }
      weight_plus    = weight_plus/double(nPdfs-1)*2.0;
      weight_minus   = weight_minus/double(nPdfs-1)*2.0;
      weight_average = weight_average/double(nPdfs);
    }

    SetValue("PdfWeightPlus_"+pdfName, weight_plus);
    SetValue("PdfWeightMinus_"+pdfName, weight_minus);
    SetValue("PdfWeightAverage_"+pdfName, weight_average);


  }
  /*
  // from example in wiki
  edm::InputTag pdfWeightTag("pdfWeights:cteq66"); // or any other PDF set
  edm::Handle<std::vector<double> > weightHandle;
  event.getByLabel(pdfWeightTag, weightHandle);
  
  std::vector<double> weights = (*weightHandle);
  std::cout << "Event weight for central PDF:" << weights[0] << std::endl;
  unsigned int nmembers = weights.size();
  for (unsigned int j=1; j<nmembers; j+=2) {
    std::cout << "Event weight for PDF variation +" << (j+1)/2 << ": " << weights[j] << std::endl;
    std::cout << "Event weight for PDF variation -" << (j+1)/2 << ": " << weights[j+1] << std::endl;
  }
  */


    
  return 0;
}
