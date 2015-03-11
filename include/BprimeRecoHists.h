#pragma once 

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"

#include "TH1F.h"

#include <string>
#include <vector>

class BprimeRecoHists: public uhh2::Hists {
 public:
  BprimeRecoHists(uhh2::Context & ctx, const std::string & dirname);
  virtual ~BprimeRecoHists();
  virtual void fill(const uhh2::Event & ev) override;
 protected:
  struct BaseHists{
    TH1F* pt, *eta, *phi, *mass; 
  };
  BaseHists book_BaseHists(const std::string & name, const std::string & label, double minPt=0, double maxPt=2000);
  template<typename T>
    void fill_BaseHists(const T & particle, BaseHists & hists, double weight);
 private: 
  //all hypothesis
  BaseHists wHad_all, wLep_all, topLep_all, topHad_all;
  //best hypothesis
  BaseHists wHad_best, wLep_best, topLep_best, topHad_best;
  //DeltaR & DeltaPhi. top referece to top -> w
  TH1F* deltaR_w_all, *deltaPhi_w_all, *deltaR_top_all;
  TH1F* deltaR_w_best, *deltaPhi_w_best, *deltaR_top_best;
  
  uhh2::Event::Handle<std::vector<BprimeContainer>> hyps;
};
