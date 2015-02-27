#pragma once

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/LorentzVector.h" 

#include <vector>

class BprimeReco :public uhh2::AnalysisModule {
 public:
  explicit BprimeReco(uhh2::Context & ctx, const std::string & label="BprimeReco");
  virtual bool process(uhh2::Event & event) override;
  bool massReco(uhh2::Event & event);
 private:
  void comb(int N, int K);
  uhh2::Event::Handle<std::vector<BprimeContainer>> hypothesis;
  uhh2::Event::Handle<FlavorParticle> primlep;
};

//more classes for TopTag and WTag Reconstruction should follow
