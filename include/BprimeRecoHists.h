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

 private:
  //string HistoNames;    
  TH1F* whad_pt, *whad_phi, *whad_eta, *whad_mass;
  TH1F* wlep_pt, *wlep_phi, *wlep_eta, *wlep_mass;
  TH1F* deltaR_w, *deltaPhi_w;

  TH1F* whad_best_pt, *whad_best_phi, *whad_best_eta, *whad_best_mass;
  TH1F* wlep_best_pt, *wlep_best_phi, *wlep_best_eta, *wlep_best_mass;
  TH1F* deltaR_w_best, *deltaPhi_w_best;

  uhh2::Event::Handle<std::vector<BprimeContainer>> hyps;
};
