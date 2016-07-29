#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"

class TopTagScalefactor: public uhh2::AnalysisModule{
 public:
  explicit TopTagScalefactor(uhh2::Context & ctx, const std::string & hyp_name);
  virtual bool process(uhh2::Event & event) override;

 private:
  uhh2::Event::Handle<BprimeContainer> recohyp;
  uhh2::Event::Handle<double> weight_toptag;
  uhh2::Event::Handle<double> weight_toptag_up;
  uhh2::Event::Handle<double> weight_toptag_down;
};
