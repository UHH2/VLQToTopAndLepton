#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/core/include/AnalysisModule.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"


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

class HTLepSelection: public uhh2::Selection{
public:
  //enum htType{HT, HtLep};
  explicit HTLepSelection(uhh2::Context & ctx, double HTLepmin);
  virtual bool passes(const uhh2::Event & event);
  
 private:
  //htType type;
  double HTLepmin;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};
class METSelection: public uhh2::Selection{
public:
  //enum htType{HT, HtLep};
  explicit METSelection(uhh2::Context & ctx, double METmin);
  virtual bool passes(const uhh2::Event & event);
  
 private:
  //htType type;
  double METmin;
};

class TwoDCut: public uhh2::Selection {
 public:
  explicit TwoDCut(float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {}
  virtual bool passes(const uhh2::Event & event) override;
  
 private:
  float min_deltaR_, min_pTrel_;
};
  /////

class ChiSquareCut: public uhh2::Selection{
 public:
  explicit ChiSquareCut(uhh2::Context & ctx, float max_chi2, float min_chi2=-1, const std::string & hyp_name="BprimeReco", int recotyp =-1);
  virtual bool passes(const uhh2::Event & event) override;

 private:
  float min_, max_;
  int recotyp_;
  uhh2::Event::Handle<BprimeContainer> recohyp;
};

class PtRatioWTCut: public uhh2::Selection{
 public:
  explicit PtRatioWTCut(uhh2::Context & ctx, float min, float max, const std::string & hyp_name="BprimeReco");
  virtual bool passes(const uhh2::Event & event) override;
 private:
   float min_, max_;
   uhh2::Event::Handle<BprimeContainer> recohyp;
};

