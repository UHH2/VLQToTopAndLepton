#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 

#include "UHH2/VLQToTopAndLepton/include/GenHT.h"

class WJets_reweight: public uhh2::AnalysisModule {
 public:
  explicit WJets_reweight(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
 private:
  uhh2::Event::Handle<double> GenHT_val;
};
