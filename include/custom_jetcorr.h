#pragma once

#include <string>

//#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 

#include "UHH2/common/include/JetCorrections.h" 
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetIds.h"



//#include "UHH2/VLQToTopAndLepton/include/Utils.h"


class custom_jetcorr {
 public:
  custom_jetcorr(uhh2::Context& ctx);
  void copy_jets(uhh2::Event & event);
  uhh2::Event::Handle<std::vector<Jet>> get_jetHandle(int i){return jet_handles[i];}
  std::string get_namestring(int i){return combination_names[i];}
  std::vector<uhh2::Event::Handle<MET>> get_methandles(){return  met_handles;}
  void apply_corr(uhh2::Event & event);
  void clean(uhh2::Event & event);
  void set_met(uhh2::Event & event);
  
  
 private:
  std::vector<std::unique_ptr<JetResolutionSmearer>> jer_correctors;
  std::vector<std::unique_ptr<JetCorrector>> jec_correctors;
  std::vector<uhh2::Event::Handle<std::vector<Jet>>> jet_handles;
  std::vector<std::unique_ptr<JetCleaner>> jet_cleaners;
  
  std::vector<uhh2::Event::Handle<MET>> met_handles;

  std::vector<std::string> combination_names = {"jer_up","jer_down","jec_up","jec_down","jec_jer_down","jec_jer_up","jec_jer_up_down","jec_jer_down_up","nominal"};
  //first entry is jer, second is jec
  std::vector<std::vector<int>> combination_ints = {{1,0},{-1,0},{0,1},{0,-1},{-1,-1},{1,1},{-1,1},{1,-1},{0,0}};
  
};
