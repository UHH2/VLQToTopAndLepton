#pragma once

#include <string>
#include <iostream>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Utils.h" 

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"

#include "UHH2/VLQToTopAndLepton/include/RunDependendJetCorr.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeDiscriminator.h"

class AK8_cor_reco :public uhh2::AnalysisModule{
 public:
  explicit AK8_cor_reco(uhh2::Context & ctx, std::string topjetcollection, TopJetId toptagid, TopJetId wjetId);
  virtual bool process(uhh2::Event & event) override;

 private:
  std::vector<std::string> correction_names = {"jec_up","jec_down"};
  //first entry is jer, second is jec
  std::vector<int> correction_int = {1,-1};
  std::vector<std::unique_ptr<GenericTopJetCorrector>> corrector_topjet;
  //std::vector<std::unique_ptr<TopJetCorrector>> corrector_topjet;
  std::vector<std::unique_ptr<GenericSubJetCorrector>> corrector_subjet;
  std::vector<uhh2::Event::Handle<std::vector<TopJet>>> topjet_handles;
  std::vector<std::unique_ptr<BprimeReco>> reco;
  std::vector<std::unique_ptr<BprimeDiscriminator>> dis;
  uhh2::Event::Handle<std::vector<TopJet> > h_topjets;
  bool is_mc=false;
};
