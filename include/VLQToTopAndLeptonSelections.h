#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/core/include/AnalysisModule.h"

class GenParticleFilter: public uhh2::Selection {
 public:
  explicit GenParticleFilter(int pdgId, int nmin, int nmax=-1);
  virtual bool passes(const uhh2::Event & event);
 private:
  int pdgId, nmin, nmax;
  

};

class HTSelection: public uhh2::Selection{
 public:
  //enum htType{HT, HtLep};
  explicit HTSelection(uhh2::Context & ctx, double HTmin);
  virtual bool passes(const uhh2::Event & event);
  
 private:
  //htType type;
  double HTmin;
  uhh2::Event::Handle<double> ht;
};

class STSelection: public uhh2::Selection{
 public:
  //enum htType{HT, HtLep};
  explicit STSelection(uhh2::Context & ctx, double STmin);
  virtual bool passes(const uhh2::Event & event);
  
 private:
  //htType type;
  double STmin;
  uhh2::Event::Handle<double> ht;
  uhh2::Event::Handle<FlavorParticle> h_primlep;

};

