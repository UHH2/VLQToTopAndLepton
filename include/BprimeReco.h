#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/LorentzVector.h"
 
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"
#include "UHH2/VLQToTopAndLepton/include/NeutrinoFit.h"

#include <vector>

class BprimeReco :public uhh2::AnalysisModule {
 public:
  explicit BprimeReco(uhh2::Context & ctx, const std::string & label="BprimeReco");
  virtual bool process(uhh2::Event & event) override;
  bool massReco(uhh2::Event & event);
  bool TopJetReco(uhh2::Event & event, double dRmin = 0);
  bool BTagReco(uhh2::Event & event);
  bool hadronicW(uhh2::Event & event, double dRmin = 0);
  void set_jetRecoId(const JetId & my_jetId){jetId= my_jetId;}
  void set_topjetRecoId(const TopJetId & my_topjetId){topjetId =my_topjetId;}
  void set_wjetRecoId(const TopJetId & my_wjetId){wjetId =my_wjetId;}
  bool set_topjetCollection(uhh2::Context & ctx, const std::string & topjetCollectionName);
  bool set_jetCollection(uhh2::Context & ctx, const std::string & jetCollectionName);  
  void set_forwardId(const JetId & forwardId){forwardJetId = forwardId;}
  void set_bjetId(const JetId & bjetId_){bjetId=bjetId_;}

 private:
  void comb(int N, int K);
  int count_bjets(std::vector<Jet> & jets, uhh2::Event & event);
  double forwardeta(std::vector<Jet> & jets, uhh2::Event & event);
  std::vector<BprimeContainer> reconstruct_WHyps(const std::vector<Jet> & jets, const std::vector<LorentzVector> & Wleps, double cutoff_WHad_min=50, double cutoff_WHad_max=400);
  template<typename T>
    bool passes_id(const T & object, const uhh2::Event & event, const boost::optional<std::function<bool (const T &, const uhh2::Event & )>> & object_id);
  uhh2::Event::Handle<std::vector<BprimeContainer>> hypothesis;
  uhh2::Event::Handle<FlavorParticle> primlep;
  uhh2::Event::Handle<std::vector<Jet>> jet_collection;
  uhh2::Event::Handle<std::vector<TopJet>> topjet_collection;
  bool topjetCollBool, jetCollBool;
  boost::optional<JetId> jetId;
  boost::optional<TopJetId> topjetId;
  boost::optional<TopJetId> wjetId;
  boost::optional<JetId> forwardJetId;
  boost::optional<JetId> bjetId;

  NeutrinoFit FitNeutrino;

};


//more classes for TopTag and WTag Reconstruction should follow
