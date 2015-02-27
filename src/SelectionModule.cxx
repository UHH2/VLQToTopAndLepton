#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/CleaningModules.h"

#include "UHH2/common/include/EventVariables.h"

#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"

#include "UHH2/common/include/TriggerSelection.h" 
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetCorrections.h" 
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TTbarReconstruction.h"




#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeRecoHists.h"


using namespace std;
using namespace uhh2;


class SelectionModule: public AnalysisModule {
public:

  explicit SelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;

  std::unique_ptr<BprimeRecoHists> RecoHists;
  std::unique_ptr<AnalysisModule> ht, lepton;
  std::unique_ptr<BprimeReco> Reco;

};



SelectionModule::SelectionModule(Context& ctx){
  Reco.reset(new BprimeReco(ctx)); 
  ht.reset(new HTCalc(ctx));
  lepton.reset(new PrimaryLepton(ctx));
  RecoHists.reset(new BprimeRecoHists(ctx,"RecoBprime"));
}

bool SelectionModule::process(Event & event){

  ht->process(event);
  lepton->process(event);
  

  if(Reco->massReco(event))
    RecoHists->fill(event);
 

  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SelectionModule)
