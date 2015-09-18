#pragma once 

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeGenContainer.h"

#include "TH1F.h"
#include "TH2F.h"

#include <string>
#include <vector>

class BprimeHypHists: public uhh2::Hists {
 public:
  BprimeHypHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyp_name);
  virtual ~BprimeHypHists();
  virtual void fill(const uhh2::Event & ev) override;
 protected:
  struct BaseHists{
    TH1F* pt, *eta, *phi, *mass; 
  };
  BaseHists book_BaseHists(const std::string & name, const std::string & label, double minMass=0, double maxMass=600, double minPt=0, double maxPt=2000);
  template<typename T>
    void fill_BaseHists(const T & particle, BaseHists & hists, double weight);
 private: 
  //hypothesis
  BaseHists wHad, wLep, topLep, topHad;
  BaseHists mass, mass_lep, mass_had;
  BaseHists ttbar;
  BaseHists forward, balance, combination, no_forwardJet;
  //DeltaR & DeltaPhi. top referece to top -> w
  TH1F* deltaR_w, *deltaPhi_w, *deltaR_top, *deltaR_wtop;
  TH1F* pTratio_wtop, *pTratio_toptop, *pTratio_ww;
  TH1F* recotype_h, *chiDis, *chiDis_lep, *chiDis_had;
  
  TH1F* wHad_res_pt, *wHad_res_E, *wHad_res_mass, *wHad_res_phi, *wHad_res_eta, *wHad_res_deltaR;
  TH1F* wLep_res_pt, *wLep_res_E, *wLep_res_mass, *wLep_res_phi, *wLep_res_eta, *wLep_res_deltaR;
  TH1F* topHad_res_pt, *topHad_res_E, *topHad_res_mass, *topHad_res_phi, *topHad_res_eta, *topHad_res_deltaR;
  TH1F* topLep_res_pt, *topLep_res_E, *topLep_res_mass, *topLep_res_phi, *topLep_res_eta, *topLep_res_deltaR;
  TH1F* Bprime_res_pt, *Bprime_res_E, *Bprime_res_mass, *Bprime_res_phi, *Bprime_res_eta, *Bprime_res_deltaR;

  TH2F* topReco_dR_pT_lep, *topReco_dR_pTres_lep, *topReco_dR_pT_had, *topReco_dR_pTres_had, *wReco_dR_pT_lep, *wReco_dR_pTres_lep, *wReco_dR_pT_had, *wReco_dR_pTres_had;

  TH2F* chi_top_pT, *chi_wlep_pT, *chi_whad_pT, chi_ST;
  TH2F* chi_deltaR_w_top;
  TH1F* deltaR_forward_B; 
  TH2F* forward_pt_eta;


  BaseHists matched_top_lep, matched_top_had, matched_W_lep, matched_W_had;
  BaseHists matched_tops, matched_top_matched_W, matched_W_matched_top, matched_top_unmatched_W, matched_W_unmatched_top;
  BaseHists semi_matched_tops; 
  
  uhh2::Event::Handle<BprimeContainer> recohyp;
  uhh2::Event::Handle<BprimeGenContainer> gen;
};
