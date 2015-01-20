#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CleaningModules.h"

#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/JetHists.h"

#include "UHH2/common/include/TriggerSelection.h" 

#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"
#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonHists.h"
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"

using namespace std;
using namespace uhh2;

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class VLQToTopAndLeptonModule: public AnalysisModule {
public:
  
  explicit VLQToTopAndLeptonModule(Context & ctx);
  virtual bool process(Event & event) override;
  
private:
  
  std::unique_ptr<JetCleaner> jetcleaner;
   
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel, bsel, trigger;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts, h_njet, h_bsel, h_ele, h_ele_trigger, h_vlqGen;

};


VLQToTopAndLeptonModule::VLQToTopAndLeptonModule(Context & ctx){
  // In the constructor, the typical tasks are to create
  // other modules like cleaners (1), selections (2) and Hist classes (3).
  // But you can do more and e.g. access the configuration, as shown below.
  
    
  // If needed, access the configuration of the module here, e.g.:
  
  
  // If running in SFrame, the keys "dataset_version", "dataset_type" and "dataset_lumi"
  // are set to the according values in the xml file. For CMSSW, these are
  // not set automatically, but can be set in the python config file.

  string testvalue = ctx.get("SelectionTriggerName", "<not set>");
  cout << "Trigger in the configuration was: " << testvalue << endl;
  
  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }
  
  //setup trigger
  trigger.reset(new TriggerSelection(testvalue));


  // 1. setup other modules. Here, only the jet cleaner
  jetcleaner.reset(new JetCleaner(30.0, 2.4));
  
  // 2. set up selections:
  njet_sel.reset(new NJetSelection(2));
  bsel.reset(new NBTagSelection(1));
  
  // 3. Set up Hists classes:
  h_vlqGen.reset(new VLQGenHists(ctx,"VLQGenHists"));
  
  h_nocuts.reset(new VLQToTopAndLeptonHists(ctx, "NoCuts"));
  h_njet.reset(new VLQToTopAndLeptonHists(ctx, "Njet"));
  h_bsel.reset(new VLQToTopAndLeptonHists(ctx, "Bsel"));
  h_ele.reset(new ElectronHists(ctx, "ele_nocuts"));
  h_ele_trigger.reset(new ElectronHists(ctx, "ele_nocuts"));
}


bool VLQToTopAndLeptonModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
    //cout << "VLQToTopAndLeptonModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    
    // 1. run all modules; here: only jet cleaning.
    jetcleaner->process(event);
    
    // 2. test selections and fill histograms
    
    h_vlqGen->fill(event);

    h_nocuts->fill(event);
    
    bool njet_selection = njet_sel->passes(event);
    if(njet_selection){
        h_njet->fill(event);
    }
    bool bjet_selection = bsel->passes(event);
    if(bjet_selection){
        h_bsel->fill(event);
    }
    h_ele->fill(event);


    if(trigger->passes(event)){
      h_ele_trigger->fill(event);

    }


    
    // 3. decide whether or not to keep the current event in the output:
    return njet_selection && bjet_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the VLQToTopAndLeptonModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(VLQToTopAndLeptonModule)
