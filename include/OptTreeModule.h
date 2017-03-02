#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/CommonModules.h"


#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"


class OptTreeModule : public uhh2::AnalysisModule {
public:
explicit OptTreeModule(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
 private:
  int numChi2, numTopTag, numTTbar;
  //Input Handles
  uhh2::Event::Handle<BprimeContainer> chi2HypHandle;
  uhh2::Event::Handle<BprimeContainer> TopTagHypHandle;
  uhh2::Event::Handle<BprimeContainer> ttbarHandle;
  uhh2::Event::Handle<BprimeContainer> wTagHypHandle;
  //Output
  uhh2::Event::Handle<double> ST;
  uhh2::Event::Handle<double> MET;
  uhh2::Event::Handle<double> HTLep;
  uhh2::Event::Handle<double> leadingLepPt;
  uhh2::Event::Handle<double> leadingJetPt;
  uhh2::Event::Handle<double> subleadJetPt;
  uhh2::Event::Handle<double> leadingTopJetPt;
  uhh2::Event::Handle<double> ForwardJetEta;
  uhh2::Event::Handle<double> numberJet;
  //uhh2::Event::Handle<double> additonalJets;
  uhh2::Event::Handle<double> chi2Mass;	
  uhh2::Event::Handle<double> chi2Val;
  uhh2::Event::Handle<double> TopTagMass;
  uhh2::Event::Handle<double> TopTagChi2;
  uhh2::Event::Handle<double> ttbarMass;
  uhh2::Event::Handle<double> ttbarChi2;
  //uhh2::Event::Handle<double> wTagMass;
  uhh2::Event::Handle<double> weight;
  uhh2::Event::Handle<double> pTW;	
  uhh2::Event::Handle<double> pTT;
  uhh2::Event::Handle<double> NBmediumtags;
  uhh2::Event::Handle<double> NBtighttags;
  uhh2::Event::Handle<double> NBloosetags;
 
//uhh2::Event::Handle<double> ;
  JetId mediumbtag_id, tightbtag_id, loosebtag_id;
};

/**
 ** To Do List for the optimizer! See above
 ** Do not forget to check weather pointer is Null
 ** Do not forget: 
 ** -> allow_null in GenericEvent
 ** -> assert(addr != nullptr) in root-utils.cxx
 */

