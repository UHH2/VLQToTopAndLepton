#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TriggerSelection.h" 
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/VLQToTopAndLepton/include/GenSelection.h"
#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"
#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeRecoHists.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeDiscriminator.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeHypHists.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeGen.h"
#include "UHH2/VLQToTopAndLepton/include/Utils.h"
#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"
#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonHists.h"
#include "UHH2/VLQToTopAndLepton/include/OptTreeModule.h"


using namespace std;
using namespace uhh2;



class OptSelModule: public AnalysisModule {
public:

  explicit OptSelModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  string Version;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<BprimeDiscriminator> chi2_combo, ttbar, cmstoptagDis, heptoptagDis;
  std::unique_ptr<AnalysisModule> ht, lepton, OptTree;
  std::unique_ptr<BprimeReco> Reco, CMSTopTagReco, HEPTopTagReco;
  std::unique_ptr<HistFactory> CutPlots;
  std::unique_ptr<MCPileupReweight> pileup_weights;
  JetId btag_medium, eta_cut, twojet, onejet, lowjet;
  MuonId muid_cut;
  TopJetId cmstopjetid, heptopjetid;
};

OptSelModule::OptSelModule(Context& ctx){
  common.reset(new CommonModules());
  common->set_jet_id(PtEtaCut(40.0,3));
  common->set_electron_id(AndId<Electron>(ElectronID_Spring15_25ns_tight_noIso, PtEtaCut(40.0, 2.1)));
  common->set_muon_id(AndId<Muon>(MuonIDTight(),PtEtaCut(40.0, 2.1)));
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->init(ctx);
  OptTree.reset(new OptTreeModule(ctx));
  cmstopjetid = AndId<TopJet>(CMSTopTag(),Tau32());
  heptopjetid = HEPTopTagV2();
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  //btag_medium = CSVBTag(CSVBTag::WP_TIGHT);
  Reco.reset(new BprimeReco(ctx));
  CMSTopTagReco.reset(new BprimeReco(ctx,"CMSTopTagReco"));
  CMSTopTagReco->set_topjetRecoId(cmstopjetid);
  HEPTopTagReco.reset(new BprimeReco(ctx,"HEPTopTagReco"));
  HEPTopTagReco->set_topjetRecoId(heptopjetid);
  HEPTopTagReco->set_topjetCollection(ctx,"patJetsHepTopTagCHSPacked_daughters");
  //Event::Handle<std::vector<TopJet>> heptopjets_handle = ctx.get_handle<std::vector<TopJet>>("patJetsHepTopTagCHSPacked_daughters");
  eta_cut = PtEtaCut(20.0,2.6);
  ht.reset(new HTCalc(ctx,eta_cut));
  lepton.reset(new PrimaryLepton(ctx));
  chi2_combo.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::chi2_combo,"BprimeReco","Chi2Dis"));
  chi2_combo->set_emptyHyp(true);
  ttbar.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::ttbar,"BprimeReco","TTbarDis"));
  ttbar->set_emptyHyp(true);
  cmstoptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::cmsTopTag,"CMSTopTagReco","CMSTopTagDis"));
  cmstoptagDis->set_emptyHyp(true);
  heptoptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::hepTopTag,"HEPTopTagReco","HEPTopTagDis"));
  heptoptagDis->set_emptyHyp(true);

  lowjet = PtEtaCut(50.0,2.4); 
  onejet = PtEtaCut(100.0, 2.4); twojet = PtEtaCut(50.0, 2.4);
  muid_cut = AndId<Muon>(MuonIDTight(), PtEtaCut(50.0, 2.4));

  CutPlots.reset(new HistFactory(ctx));
  CutPlots->setEffiHistName("CutsOpt");
  CutPlots->addSelection(make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*"),"muonTrigger");
  CutPlots->addSelection(make_unique<NElectronSelection>(0,0),"0_eleCut");
  CutPlots->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  CutPlots->addSelection(make_unique<NJetSelection>(1,-1,lowjet),"50GeV_JetCut");
  CutPlots->addSelection(make_unique<TwoDCut>(0.4,25),"2DCut");
  //CutPlots->addSelection(make_unique<ForwardJetPtEtaCut>(1.5,40),"ForwardJetCut");
  CutPlots->addSelection(make_unique<HTLepSelection>(ctx,100),"100_HTLep");
  CutPlots->addSelection(make_unique<NJetSelection>(1,-1,onejet),"100GeV_JetCut");
  
  CutPlots->addHists("ElectronHists","ElectronHists");
  CutPlots->addHists("MuonHists","MuonHists");
  CutPlots->addHists("EventHists","EventHists");
  CutPlots->addHists("JetHists","JetHists");
  CutPlots->addHists("TopJetHists","TopJetHists");
  CutPlots->addHists("VLQGenHists","VLQGenHists");
  CutPlots->addHists("TopTag_JetHists",cmstopjetid);
  CutPlots->addHists("BTagged_JetHists",btag_medium);
}

bool OptSelModule::process(Event & event){
  if(!common->process(event)) return false;
  ht->process(event);
  lepton->process(event);
  if(!CutPlots->passAndFill(event))return false;
  Reco->massReco(event);
  ttbar->process(event);
  chi2_combo->process(event);
  HEPTopTagReco->TopJetReco(event,2);
  heptoptagDis->process(event);
  CMSTopTagReco->TopJetReco(event,2);
  cmstoptagDis->process(event);
  OptTree->process(event);
  return true;
}
// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(OptSelModule)
