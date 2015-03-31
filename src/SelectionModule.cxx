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

#include "UHH2/VLQToTopAndLepton/include/GenSelection.h"
#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"
#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeRecoHists.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeGen.h"

#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"


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
  std::unique_ptr<BprimeRecoHists> Reco_wHad, Reco_wLep;
  std::unique_ptr<BprimeGen> Gen;
  std::unique_ptr<AnalysisModule> jetlepclean; 
  std::unique_ptr<HistFactory> stdPlots, TopTagPlots, BTagPlots;
  std::unique_ptr<HistFactory> TopTagBTag, TopTagNoBTag, NoTopTagBTag, NoTopTagNoBtag;
  std::unique_ptr<uhh2::Selection> twoDcut;
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<Selection> topLep,topHad;
  JetId btag_medium;
  TopJetId topjetid;
};

SelectionModule::SelectionModule(Context& ctx){
  topjetid = CMSTopTag();
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  jetlepclean.reset(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));
  Reco.reset(new BprimeReco(ctx));
  Gen.reset(new BprimeGen(ctx)); 
  ht.reset(new HTCalc(ctx));
  lepton.reset(new PrimaryLepton(ctx));
  twoDcut.reset(new TwoDCut(.4, 25.));
  jetcleaner.reset(new JetCleaner(50.0, 2.4));
 
  RecoHists.reset(new BprimeRecoHists(ctx,"RecoBprime"));

  stdPlots.reset(new HistFactory(ctx));
  stdPlots->setEffiHistName("2D");
  stdPlots->addSelection(make_unique<TwoDCut>(0.4,25),"2DCut");
  stdPlots->addHists("ElectronHists","std_ElectronHists");
  stdPlots->addHists("MuonHists","std_MuonHists");
  stdPlots->addHists("EventHists","std_EventHists");
  stdPlots->addHists("JetHists","std_JetHists");
  stdPlots->addHists("TopJetHists","std_TopJetHists");
  stdPlots->addHists("VLQGenHists","std_VLQGenHists");

  TopTagPlots.reset(new HistFactory(ctx));
  TopTagPlots->setEffiHistName("TopTag");
  TopTagPlots->addSelection(make_unique<NTopJetSelection>(1,-1,topjetid),"TopTag");
  TopTagPlots->addSelection(make_unique<NTopJetSelection>(0,0,topjetid),"AntiTopTag");
  TopTagPlots->addHists("ElectronHists","TopTag_ElectronHists");
  TopTagPlots->addHists("MuonHists","TopTag_MuonHists");
  TopTagPlots->addHists("EventHists","TopTag_EventHists");
  TopTagPlots->addHists("JetHists","TopTag_JetHists");
  TopTagPlots->addHists("TopJetHists","TopTag_TopJetHists");
  TopTagPlots->addHists("VLQGenHists","TopTag_VLQGenHists");
  TopTagPlots->addHists("BprimeRecoHists","TopTag_BprimeRecoHists");

  BTagPlots.reset(new HistFactory(ctx));
  BTagPlots->setEffiHistName("BTag");
  BTagPlots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"BTag");
  BTagPlots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"2_BTags");
  BTagPlots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"AntiBTag");
  BTagPlots->addHists("ElectronHists","BTag_ElectronHists");
  BTagPlots->addHists("MuonHists","BTag_MuonHists");
  BTagPlots->addHists("EventHists","BTag_EventHists");
  BTagPlots->addHists("JetHists","BTag_JetHists");
  BTagPlots->addHists("TopJetHists","BTag_TopJetHists");
  BTagPlots->addHists("VLQGenHists","BTag_VLQGenHists");
  BTagPlots->addHists("BprimeRecoHists","BTag_BprimeRecoHists");

  vector<int> topLepIds {6,24,13};
  vector<int> topHadIds {6,24,-54321};
  topLep.reset(new GenFamilySelection(topLepIds,2));
  topHad.reset(new GenFamilySelection(topHadIds,2));;
  Reco_wHad.reset(new BprimeRecoHists(ctx, "Gen_wHad"));
  Reco_wLep.reset(new BprimeRecoHists(ctx, "Gen_wLep"));
}

bool SelectionModule::process(Event & event){
  jetlepclean->process(event);
  ht->process(event);
  jetcleaner->process(event);
  lepton->process(event);
  stdPlots->passAndFill(event);

  Gen->process(event);
  if(!Reco->massReco(event) || !twoDcut->passes(event)) return false;
 
  if(topLep->passes(event))  Reco_wHad->fill(event);
  if(topHad->passes(event))  Reco_wLep->fill(event);

  RecoHists->fill(event);
  TopTagPlots->passAndFill(event,1);
  BTagPlots->passAndFill(event,1);

  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SelectionModule)
