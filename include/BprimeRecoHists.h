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

class BprimeRecoHists: public uhh2::Hists {
 public:
  BprimeRecoHists(uhh2::Context & ctx, const std::string & dirname);
  virtual ~BprimeRecoHists();
  virtual void fill(const uhh2::Event & ev) override;
 protected:
  struct BaseHists{
    TH1F* pt, *eta, *phi, *mass; 
  };
  BaseHists book_BaseHists(const std::string & name, const std::string & label, double minMass=0, double maxMass=600, double minPt=0, double maxPt=2000);
  template<typename T>
    void fill_BaseHists(const T & particle, BaseHists & hists, double weight);
 private: 
  //double calc_Chi(LorentzVector top, LorentzVector w);
  //all hypothesis
  BaseHists wHad_all, wLep_all, topLep_all, topHad_all;
  //best hypothesis
  BaseHists wHad_best, wLep_best, topLep_best, topHad_best;
  //BaseHists wHad_best_lepHyp, wLep_best_lepHyp, topLep_best_lepHyp, topHad_best_lepHyp;
  //BaseHists wHad_best_hadHyp, wLep_best_hadHyp, topLep_best_hadHyp, topHad_best_hadHyp;
  BaseHists Bprime, Bprime_lep, Bprime_had;	
  //DeltaR & DeltaPhi. top referece to top -> w
  TH1F* deltaR_w_all, *deltaPhi_w_all, *deltaR_top_all;
  TH1F* deltaR_w_best, *deltaPhi_w_best, *deltaR_top_best;
  TH1F* recotype_h, *chiDis, *chiDis_lep, *chiDis_had;
  TH1F* ttbar_chi;

  TH1F* wHad_res_pt, *wHad_res_E, *wHad_res_mass, *wHad_res_phi, *wHad_res_eta, *wHad_res_deltaR;
  TH1F* wLep_res_pt, *wLep_res_E, *wLep_res_mass, *wLep_res_phi, *wLep_res_eta, *wLep_res_deltaR;
  TH1F* topHad_res_pt, *topHad_res_E, *topHad_res_mass, *topHad_res_phi, *topHad_res_eta, *topHad_res_deltaR;
  TH1F* topLep_res_pt, *topLep_res_E, *topLep_res_mass, *topLep_res_phi, *topLep_res_eta, *topLep_res_deltaR;
  
  TH2F* topReco_dR_pT_lep, *topReco_dR_pTres_lep, *topReco_dR_pT_had, *topReco_dR_pTres_had, *wReco_dR_pT_lep, *wReco_dR_pTres_lep, *wReco_dR_pT_had, *wReco_dR_pTres_had;
  

  uhh2::Event::Handle<std::vector<BprimeContainer>> hyps;
  uhh2::Event::Handle<BprimeGenContainer> gen;
};
