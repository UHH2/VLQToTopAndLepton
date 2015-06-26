#pragma once

#include "UHH2/VLQToTopAndLepton/include/BprimeContainer.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/LorentzVector.h" 
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"

#include <vector>

class BprimeReco :public uhh2::AnalysisModule {
 public:
  explicit BprimeReco(uhh2::Context & ctx, const std::string & label="BprimeReco");
  virtual bool process(uhh2::Event & event) override;
  bool massReco(uhh2::Event & event);
  bool TopJetReco(uhh2::Event & event, double dRmin = 0);
  void set_jetRecoId(boost::optional<JetId> my_jetId){jetId= my_jetId;}
  void set_topjetRecoId(boost::optional<TopJetId> my_topjetId){topjetId =my_topjetId;}
  bool set_topjetCollection(uhh2::Context & ctx, const std::string & topjetCollectionName);
  bool set_jetCollection(uhh2::Context & ctx, const std::string & jetCollectionName);  
 private:
  void comb(int N, int K);

  template<typename T>
    bool passes_id(const T & object, const uhh2::Event & event, const boost::optional<std::function<bool (const T &, const uhh2::Event & )>> & object_id);
  uhh2::Event::Handle<std::vector<BprimeContainer>> hypothesis;
  uhh2::Event::Handle<FlavorParticle> primlep;
  uhh2::Event::Handle<std::vector<Jet>> jet_collection;
  uhh2::Event::Handle<std::vector<TopJet>> topjet_collection;
  bool topjetCollBool, jetCollBool;
  boost::optional<JetId> jetId;
  boost::optional<TopJetId> topjetId;
};

//more classes for TopTag and WTag Reconstruction should follow
