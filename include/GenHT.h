#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"


class GenHT: public uhh2::AnalysisModule {
 public:
  explicit GenHT(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
 private:
  uhh2::Event::Handle<double> GenHT_val;
};
