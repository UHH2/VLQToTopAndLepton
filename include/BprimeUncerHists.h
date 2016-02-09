#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeGenContainer.h"

#include "TH1F.h"

#include <vector>
#include <string>

class BprimeUncerHists: public uhh2::Hists {
 public:
  BprimeUncerHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyp_name);
  virtual ~BprimeUncerHists();
  virtual void fill(const uhh2::Event & ev) override;
 private:
  uhh2::Event::Handle<BprimeContainer> recohyp;
  std::vector<TH1F*> scale;
  std::vector<int> scale_entries = {1,2,3,4,6,8};
  std::vector<TH1F*> pdf;
};
