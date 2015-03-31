#pragma once

#include "UHH2/VLQToTopAndLepton/include/BprimeGenContainer.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/LorentzVector.h" 

#include <vector>

class BprimeGen :public uhh2::AnalysisModule {
 public:
  explicit BprimeGen(uhh2::Context & ctx, const std::string & label="BprimeGen");
  virtual bool process(uhh2::Event & event) override;
 private:
  uhh2::Event::Handle<BprimeGenContainer> BprimeGenLevel;
  bool family(std::vector<int> ties, const GenParticle & part, uhh2::Event & event);
};


