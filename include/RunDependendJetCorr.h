#pragma once

//#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/common/include/JetCorrections.h"

#include <string>

class RunDependendJetCorr :public uhh2::AnalysisModule {
 public:
  
  explicit RunDependendJetCorr(uhh2::Context & ctx, std::string topjetcollection, std::string subjetcollection="", int direction=0);
  virtual bool process(uhh2::Event & event) override;

 private:
  std::unique_ptr<GenericTopJetCorrector> jet_corrector_MC, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
  std::unique_ptr<SubJetCorrector> subjet_corrector_MC, subjet_corrector_BCD, subjet_corrector_EFearly, subjet_corrector_FlateG, subjet_corrector_H;
  std::unique_ptr<GenericSubJetCorrector> subjet_corrector_MC_coll;
  
  const int runnr_BCD = 276811;
  const int runnr_EFearly = 278802;
  const int runnr_FlateG = 280385;
  bool is_mc, subjetcollbool;
};

namespace JERFiles {

  extern const std::vector<std::string> Summer16_23Sep2016_V4_L123_AK8PFpuppi_MC;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK8PFpuppi_DATA;  
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK8PFpuppi_DATA; 
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK8PFpuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK8PFpuppi_DATA;
}
