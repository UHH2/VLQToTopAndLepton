#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeGenContainer.h"

#include <vector>
#include <string>
#include <iostream>

class BprimeDiscriminator :public uhh2::AnalysisModule {
 public:
  enum discriminatorType{ttbar, chi2_combo, lepTop, hadTop, cmsTopTag,hepTopTag};
  explicit BprimeDiscriminator(uhh2::Context & ctx, discriminatorType dis_, const std::string& RecoLabel="", const std::string Outputname="" ,const std::string& GenLabel="");
  virtual bool process(uhh2::Event & event) override;
 private:
  discriminatorType dis;
  uhh2::Event::Handle<std::vector<BprimeContainer>> hyps;
  uhh2::Event::Handle<BprimeGenContainer> gen;
  BprimeContainer ttbar_dis(uhh2::Event & event);
  BprimeContainer chiCombo_dis(uhh2::Event & event);
  BprimeContainer cmsTopTag_dis(uhh2::Event & event);
  uhh2::Event::Handle<BprimeContainer> resultHyp;
};
