#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"


class OptTreeModule : public uhh2::AnalysisModule {
public:
explicit OptTreeModule(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
 private:
  uhh2::Event::Handle<double> leadingJetPt;
  uhh2::Event::Handle<double> numberJet;
/*
  uhh2::Event::Handle<double> chi2Mass;	
  uhh2::Event::Handle<double> chi2Val;
  uhh2::Event::Handle<double> pTW;	
  uhh2::Event::Handle<double> pTT;    
  uhh2::Event::Handle<double> ttbarChi2;
  uhh2::Event::Handle<double> cmsTopTagMass;
  uhh2::Event::Handle<double> cmsTopTagChi2;
*/ 
//uhh2::Event::Handle<double> ;
};

/**
 ** To Do List for the optimizer! See above
 ** Do not forget to check weather pointer is Null
 ** Do not forget: 
 ** -> allow_null in GenericEvent
 ** -> assert(addr != nullptr) in root-utils.cxx
 */
