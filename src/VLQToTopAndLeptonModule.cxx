#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 


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



#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"
#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonHists.h"
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"
#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"

#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"


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

  
  std::unique_ptr<HistFactory> muonFactory, eleFactory;

  // object cleaner
  std::unique_ptr<AnalysisModule> ht;
  std::unique_ptr<AnalysisModule>  jetCorr; 
  std::unique_ptr<ElectronCleaner> elecleaner;
  std::unique_ptr<MuonCleaner> mucleaner;  
  std::unique_ptr<JetCleaner> jetcleaner;

  //special ids btag, toptag...
  JetId btag_medium;


  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
  std::unique_ptr<Selection> njetSel, trigger, muonSel, eleSel;

  std::unique_ptr<Selection> muonTrigger, muonMuon, muonJets, muonTop;
  std::unique_ptr<Selection> eleTrigger, eleEle, eleJets, eleTop;

  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_vlqGen;


  //no cuts
  std::unique_ptr<Hists> h_jets_noCuts, h_muon_noCuts, h_ele_noCuts, h_event_noCuts;
  //Lepton Selection
  std::unique_ptr<Hists> h_jets_leptonCuts, h_muon_leptonCuts, h_ele_leptonCuts, h_event_leptonCuts;
  //Jets Selection
  std::unique_ptr<Hists> h_jets_jetCuts, h_muon_jetCuts, h_ele_jetCuts, h_event_jetCuts;
  //Lepton+Jet Selection
  std::unique_ptr<Hists> h_jets_jetLepCuts, h_muon_jetLepCuts, h_ele_jetLepCuts, h_event_jetLepCuts;

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
  
 
  // 1. setup other modules. Here, only the jet cleaner
  jetcleaner.reset(new JetCleaner(30.0, 2.4));
  
  // 2. set up selections:
  
  jetCorr.reset(new JetCorrector(ctx,JERFiles::PHYS14_L123_MC));
  ht.reset(new HTCalc(ctx));
  elecleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_Spring15_50ns_medium_noIso, PtEtaCut(20.0, 2.4))));
  mucleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(),PtEtaCut(20.0, 2.1))));

  muonSel.reset(new NMuonSelection(1));
  eleSel.reset(new NElectronSelection(1));
  njetSel.reset(new NJetSelection(2));

  //HLT_Mu40_v1 HLT_Ele32_eta2p1_WP85_Gsf_v1

  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);

  muonTrigger.reset(new TriggerSelection("HLT_Mu40_v*")); muonMuon.reset(new NMuonSelection(1)); muonJets.reset(new NJetSelection(2));muonTop.reset(new NTopJetSelection(1));
  eleTrigger.reset(new TriggerSelection("HLT_Ele32_eta2p1_WP85_Gsf_v*")); eleEle.reset(new NElectronSelection(1)); eleJets.reset(new NJetSelection(2));eleTop.reset(new NTopJetSelection(1));

  muonFactory.reset(new HistFactory(ctx,"muonEffis.txt"));
  muonFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"muonTrigger");
  muonFactory->addSelection(move(muonMuon),"muonChannel_muonCut");
  muonFactory->addSelection(move(muonJets),"muonChannel_JetCut");
  muonFactory->addSelection(move(muonTop), "muonChannel_TopJetCut");
  muonFactory->addSelection(make_unique<NJetSelection>(2,-1, btag_medium), "muonChannel_BTagMedium");
 
  muonFactory->addHists("ElectronHists","muonChannel_ElectronHists");
  muonFactory->addHists("MuonHists","muonChannel_MuonHists");
  muonFactory->addHists("EventHists","muonChannel_EventHists");
  muonFactory->addHists("JetHists","muonChannel_JetHists");
  muonFactory->addHists("TopJetHists","muonChannel_TopJetHists");


  eleFactory.reset(new HistFactory(ctx,"eleEffis.txt"));
  eleFactory->addSelection(move(eleTrigger),"eleChannel_Trigger");
  eleFactory->addSelection(move(eleEle),"eleChannel_eleCut");
  eleFactory->addSelection(move(eleJets),"eleChannel_JetCut");
  eleFactory->addSelection(move(eleTop),"eleChannel_TopJetCut");
  eleFactory->addSelection(make_unique<NJetSelection>(2,-1, btag_medium), "eleChannel_BTagMedium");

  eleFactory->addHists("ElectronHists","eleChannel_ElectronHists");
  eleFactory->addHists("MuonHists","eleChannel_MuonHists");
  eleFactory->addHists("EventHists","eleChannel_EventHists");
  eleFactory->addHists("JetHists","eleChannel_JetHists");
  eleFactory->addHists("TopJetHists","eleChannel_TopJetHists");
  

  // 3. Set up Hists classes:
  h_vlqGen.reset(new VLQGenHists(ctx,"VLQGenHists"));
  
  //no cuts
  h_ele_noCuts.reset(new ElectronHists(ctx, "ele_noCuts"));
  h_jets_noCuts.reset(new JetHists(ctx, "jets_noCuts"));
  h_muon_noCuts.reset(new MuonHists(ctx, "muon_noCuts"));
  h_event_noCuts.reset(new EventHists(ctx,"event_noCuts"));
  
  //Lep cut
  h_ele_leptonCuts.reset(new ElectronHists(ctx, "ele_leptonCuts"));
  h_jets_leptonCuts.reset(new JetHists(ctx, "jets_leptonCuts"));
  h_muon_leptonCuts.reset(new MuonHists(ctx, "muon_leptonCuts"));
  h_event_leptonCuts.reset(new EventHists(ctx,"event_leptonCuts"));

  //jet cut
  h_ele_jetCuts.reset(new ElectronHists(ctx, "ele_jetCuts"));
  h_jets_jetCuts.reset(new JetHists(ctx, "jets_jetCuts"));
  h_muon_jetCuts.reset(new MuonHists(ctx, "muon_jetCuts"));
  h_event_jetCuts.reset(new EventHists(ctx,"event_jetCuts"));

  //lep+jet cuts
  h_ele_jetLepCuts.reset(new ElectronHists(ctx, "ele_jetLepCuts"));
  h_jets_jetLepCuts.reset(new JetHists(ctx, "jets_jetLepCuts"));
  h_muon_jetLepCuts.reset(new MuonHists(ctx, "muon_jetLepCuts"));
  h_event_jetLepCuts.reset(new EventHists(ctx,"event_jetLepCuts"));

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

    jetCorr->process(event);
    jetcleaner->process(event);
    elecleaner->process(event);
    mucleaner->process(event);

    ht->process(event);
    //trigger->passes(event);

    // 2. test selections and fill histograms
    
    h_vlqGen->fill(event);

    h_ele_noCuts->fill(event);
    h_jets_noCuts->fill(event);
    h_muon_noCuts->fill(event);
    h_event_noCuts->fill(event);


    if(muonSel->passes(event) || eleSel->passes(event)){
      h_ele_leptonCuts->fill(event);
      h_jets_leptonCuts->fill(event);
      h_muon_leptonCuts->fill(event);
      h_event_leptonCuts->fill(event);
    }
    
    if(njetSel->passes(event)){
      h_ele_jetCuts->fill(event);
      h_jets_jetCuts->fill(event);
      h_muon_jetCuts->fill(event);
      h_event_jetCuts->fill(event);
    }

    if((muonSel->passes(event) || eleSel->passes(event)) && njetSel->passes(event)){
      h_ele_jetLepCuts->fill(event);
      h_jets_jetLepCuts->fill(event);
      h_muon_jetLepCuts->fill(event);
      h_event_jetLepCuts->fill(event);
    }
    



    muonFactory->passAndFill(event);
    eleFactory->passAndFill(event);


    
    // 3. decide whether or not to keep the current event in the output:
    return true;//njet_selection && bjet_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the VLQToTopAndLeptonModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(VLQToTopAndLeptonModule)
