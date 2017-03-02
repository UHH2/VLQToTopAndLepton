#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 

#include "TF1.h"

class WJetsReweight: public uhh2::AnalysisModule {
 public:
  explicit WJetsReweight(uhh2::Context & ctx, std::string mode_="central");
  virtual bool process(uhh2::Event & event) override;
 private:
  bool work;
  TF1 *wpt;
  std::string mode;
  uhh2::Event::Handle<float> h_wjets_reweight;
  uhh2::Event::Handle<float> h_wjets_reweight_up;
  uhh2::Event::Handle<float> h_wjets_reweight_down;
  double plin[2]={1.46292,-8.99954e-04};
  double plin_up[2]={1.46292+3.11326e-02,-8.99954e-04+1.37306e-04};
  double plin_down[2]={1.46292-3.11326e-02,-8.99954e-04-1.37306e-04};
};
