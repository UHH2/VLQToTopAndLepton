#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/CommonModules.h"
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
#include "UHH2/common/include/TTbarReconstruction.h"


#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"
#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonHists.h"
#include "UHH2/VLQToTopAndLepton/include/GenSelection.h"
//#include "UHH2/VLQToTopAndLepton/include/HTSelection.h"
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"
#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"
#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"
#include "UHH2/VLQToTopAndLepton/include/Utils.h"

using namespace std;
using namespace uhh2;


class TriggerModule: public AnalysisModule {
public:

  explicit TriggerModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> muonTrigger, eleTrigger, dataTrigger;
  std::unique_ptr<AndSelection>  topHadSel, topLepSel;

  JetId btag_medium;
  MuonId muid_cut;
  JetId twojet, onejet, lowjet;
  TopJetId topjet;
  //std::unique_ptr<BprimeReco> Reco;
};



TriggerModule::TriggerModule(Context& ctx){
  //Reco.reset(new BprimeReco(ctx)); 
  //Version  = ctx.get("dataset_version", "<not set>");
  
  common.reset(new CommonModules());
  common->set_jet_id(PtEtaCut(40.0,2.7));
  common->set_electron_id(AndId<Electron>(ElectronID_PHYS14_25ns_tight_noIso, PtEtaCut(20.0, 2.1)));
  common->set_muon_id(AndId<Muon>(MuonIDTight(),PtEtaCut(20.0, 2.1)));
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->disable_mclumiweight();
  common->disable_mcpileupreweight();
  common->init(ctx);
  lepton.reset(new PrimaryLepton(ctx));
  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));


  dataTrigger.reset(new HistFactory(ctx));
  dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*"),"HLT_Mu45_eta2p1");					     
  dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu50_v*"),"HLT_Mu50");						     
  dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_v*"),"HLT_IsoMu24_eta2p1");			     
  dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v*"),"HLT_Mu40_eta2p1_PFJet200_PFJet50");
 
  muonTrigger.reset(new HistFactory(ctx));
  //muonTrigger.reset(new HistFactory(ctx));
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*"),"HLT_Mu45_eta2p1");					     
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu50_v*"),"HLT_Mu50");						     
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_v*"),"HLT_IsoMu24_eta2p1");			     
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v*"),"HLT_Mu40_eta2p1_PFJet200_PFJet50");
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*"),make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_v*")),"HLT_OR_Mu45_IsoMu24");
  muonTrigger->addHists("MuonHists","trigger_MuonHists");
  muonTrigger->addHists("JetHists","trigger_JetHists");

  vector<int> topLep {6,24,13};
  vector<int> topHad {6,24,-54321};
  vector<int> wHad {24,-54321};
  vector<int> wLep {24,13};

  topHadSel.reset(new AndSelection(ctx));
  topLepSel.reset(new AndSelection(ctx));
  
  topHadSel->add("BToThad",make_unique<GenFamilySelection>(topHad,2));
  topHadSel->add("BToWlep",make_unique<GenFamilySelection>(wLep,2));
  topLepSel->add("BToTlep",make_unique<GenFamilySelection>(topLep,2));
  topLepSel->add("BToWhad",make_unique<GenFamilySelection>(wHad,2));
 
}

bool TriggerModule::process(Event & event){
  if(!common->process(event)) return false;
  lepton->process(event);
  if(!event.isRealData){
    //if(topHadSel->passes(event) || topLepSel->passes(event))
    vlqGenHists->fill(event);
    muonTrigger->passAndFill(event);
  }
  else{
    dataTrigger->passAndFill(event);
  }

  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TriggerModule)
