#pragma once


#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/core/include/Hists.h"

#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

#include <vector>

/**
 *   Example class for booking and filling histograms, the new version using AnalysisModule mechanisms.
 */



class Lepton2DHist: public uhh2::Hists {
 public:
  // use the same constructor arguments as Hists for forwarding:
  Lepton2DHist(uhh2::Context & ctx, const std::string & dirname);
  
  virtual void fill(const uhh2::Event & ev) override;
  virtual ~Lepton2DHist();
  
 private:
  
  TH2F* leadingPT_muon_pt_eta, *leadingPT_electron_pt_eta;
  TH2F* secondPT_muon_pt_eta, *secondPT_electron_pt_eta;

};
