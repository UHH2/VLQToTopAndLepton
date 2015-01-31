#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"

#include "UHH2/common/include/CleaningModules.h"

#include "UHH2/common/include/EventVariables.h"

#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"

#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"

#include "UHH2/common/include/TriggerSelection.h" 
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetCorrections.h" 
#include "UHH2/common/include/JetIds.h"
//#include "UHH2/common/include/HTCalc.h"


#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"
#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonHists.h"
#include "UHH2/VLQToTopAndLepton/include/GenSelection.h"

#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"

#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"


using namespace std;
using namespace uhh2;


class GenTestModule: public AnalysisModule {
public:

  explicit GenTestModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  std::unique_ptr<HistFactory> muonFactory, eleFactory;
  
  std::unique_ptr<Hists> h_vlqGen;

  std::unique_ptr<AnalysisModule>  jetCorr; 
  std::unique_ptr<ElectronCleaner> elecleaner;
  std::unique_ptr<MuonCleaner> mucleaner;  
  std::unique_ptr<JetCleaner> jetcleaner;


  std::unique_ptr<Selection> HiggsFilter, ZFilter;
  std::unique_ptr<Selection> muonSel, eleSel;

  std::unique_ptr<Selection>  familySel;


};



GenTestModule::GenTestModule(Context& ctx){

  jetcleaner.reset(new JetCleaner(30.0, 2.4));
  jetCorr.reset(new JetCorrector(JERFiles::PHYS14_L123_MC));
  elecleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_CSA14_50ns_medium, PtEtaCut(25.0, 2.4))));
  mucleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(),PtEtaCut(25.0, 2.1))));


  muonSel.reset(new NMuonSelection(1));
  eleSel.reset(new NElectronSelection(1));

  vector<int> topMuonFamily {6,24,13};
  familySel.reset(new GenFamilySelection(topMuonFamily,2));

  HiggsFilter.reset(new GenParticleFilter(25,0,0));
  ZFilter.reset(new GenParticleFilter(23,0,0));

  h_vlqGen.reset(new VLQGenHists(ctx,"VLQGenHists"));
}




bool GenTestModule::process(Event & event){

    jetCorr->process(event);
    jetcleaner->process(event);
    elecleaner->process(event);
    mucleaner->process(event);

    if(!HiggsFilter->passes(event)||!ZFilter->passes(event)) return false;

    //if(familySel->passes(event))cout<<"passed"<<endl;
    //else cout<<"fail"<<endl;
    //event.genparticles->at(0).Print(event.genparticles);

    h_vlqGen->fill(event);


    if(muonSel->passes(event)){


    }
    else if(eleSel->passes(event)){


    }
    
  

    return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(GenTestModule)
