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


class EventKinematicHists: public uhh2::Hists{
 public:
  explicit EventKinematicHists(uhh2::Context & ctx, const std::string & dirname, std::string hyp_name ="");
  virtual ~EventKinematicHists();
  virtual void fill(const uhh2::Event & ev) override;
 protected:
  struct BaseHists{
    TH1F* pt, *eta, *phi, *mass, *energy; 
  };
  BaseHists book_BaseHists(const std::string & name, const std::string & label, double minMass=0, double maxMass=600, double minPt=0, double maxPt=2000, double minE=0, double maxE=2000);
  template<typename T>
    void fill_BaseHists(const T & particle, BaseHists & hists, double weight);

 private:
  TH1F *leadingAk4_lepton_dr, *leadingAk4_lepton_dphi;
  TH1F* leadingAk8_lepton_dr, *leadingAk8_lepton_dphi, *massleadingAk8_lepton_dr, *massleadingAk8_lepton_dphi;
  TH1F* leadingEtadrmin;
  TH2F* leadingEtadrmin_eta, *leadingEtadrmin_energy, *energy_eta;
  TH1F* massleadingAk8_prunedmass;
  TH1F* centrality_ak4, *centrality_ak8;
  TH2F* leadingAk4_lepton_dphi_mass;
  TH1F* noforward_Jet, *forwardjet_iso;
  BaseHists leadingEtaJetHists, secondEtaJetHists, leadingAk8Hists, leadingAk4Hists, lbHists, METlHists, recoForward;
  uhh2::Event::Handle<FlavorParticle> primlep;
  uhh2::Event::Handle<BprimeContainer> recohyp;
  std::string hyp_name_;
};
