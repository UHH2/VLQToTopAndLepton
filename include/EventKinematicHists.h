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
  EventKinematicHists(uhh2::Context & ctx, const std::string & dirname);
  virtual ~EventKinematicHists();
  virtual void fill(const uhh2::Event & ev) override;
 protected:
  struct BaseHists{
    TH1F* pt, *eta, *phi, *mass; 
  };
  BaseHists book_BaseHists(const std::string & name, const std::string & label, double minMass=0, double maxMass=600, double minPt=0, double maxPt=2000);
  template<typename T>
    void fill_BaseHists(const T & particle, BaseHists & hists, double weight);

 private:
  TH1F* leadingAk8_lepton_dr, *leadingAk4_lepton_dr;
  TH1F* leadingAk8_lepton_dphi, *leadingAk4_lepton_dphi;
  BaseHists leadingEtaJetHists, secondEtaJetHists, leadingAk8Hists, leadingAk4Hists;
  uhh2::Event::Handle<FlavorParticle> primlep;
};
