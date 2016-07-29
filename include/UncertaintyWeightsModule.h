#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include <vector>

class UncertaintyWeightsModule: public uhh2::AnalysisModule{
 public:
  explicit UncertaintyWeightsModule(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;

 private:
  uhh2::Event::Handle<double> scaleWeight_up;
  uhh2::Event::Handle<double> scaleWeight_down;
  uhh2::Event::Handle<double> pdfWeight;

  std::vector<int> scale_entries = {1,2,3,4,6,8};
};




