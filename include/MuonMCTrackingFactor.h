#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "TGraphAsymmErrors.h"


#include <string>

class MuonMCTrackingFactor: public AnalysisModule {
 public:
  MuonMCTrackingFactor(std::string filedir_, std::string graph_name_, int type_, std::string name="muon_track");
  virtual bool process(Event & event) override;
 private:

  TGraphAsymmErrors* eff;
  
  uhh2::Event::Handle<double> weight_nom;
  uhh2::Event::Handle<double> weight_up;
  uhh2::Event::Handle<double> weight_down;
};

   
