#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 

#include "UHH2/VLQToTopAndLepton/include/GenHT.h"

#include "TF1.h"

class WJetsReweight: public uhh2::AnalysisModule {
 public:
  explicit WJetsReweight(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
 private:
  bool work;
  TF1* wjets_inc, *wjets_ht, *wpt;
};
