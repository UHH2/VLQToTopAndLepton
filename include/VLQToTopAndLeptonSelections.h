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
  explicit METSelection(double METmin);
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

class PTWhadCut: public uhh2::Selection{
 public:
  explicit PTWhadCut(uhh2::Context & ctx, float min, float max=-1, const std::string & hyp_name="BprimeReco");
  virtual bool passes(const uhh2::Event & event) override;
 private:
   float min_, max_;
   uhh2::Event::Handle<BprimeContainer> recohyp;
};

class ForwardJetPtEtaCut: public uhh2::Selection{
 public:
  explicit ForwardJetPtEtaCut(uhh2::Context & ctx, float minEta, float maxEta=-1, float minPt=0, float maxPt=-1, float minDrmin = -1, float minEnergy =-1, std::string hyp_name ="");
  virtual bool passes(const uhh2::Event & event) override;
 private:
  float minEta_, maxEta_, minPt_, maxPt_, minDrmin_, minEnergy_;
  uhh2::Event::Handle<BprimeContainer> recohyp;
  std::string hyp_name_;
};


class NSubJetCut: public uhh2::Selection{
 public:
  explicit NSubJetCut(int min_subjets_, int max_subjets_=2, int first_topjet_ =1, int last_topjet_=1);
  virtual bool passes(const uhh2::Event & event) override;
 private:
  int min_subjets, max_subjets, first_topjet, last_topjet;
};

class TopJetMassCut: public uhh2::Selection{
 public:
  explicit TopJetMassCut(float mass_);
  virtual bool passes(const uhh2::Event & event) override;
 private:
  float mass;
};
