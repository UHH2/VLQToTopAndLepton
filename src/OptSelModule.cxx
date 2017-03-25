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
#include "UHH2/VLQToTopAndLepton/include/WTopJet.h"

using namespace std;
using namespace uhh2;



class OptSelModule: public AnalysisModule {
public:

  explicit OptSelModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  string Version;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<BprimeDiscriminator> chi2_combo, ttbar, toptagDis, heptoptagDis,chi2_wtag;
  std::unique_ptr<AnalysisModule>  lepton, OptTree;
  std::unique_ptr<BprimeReco> Reco, TopTagReco, HEPTopTagReco, WTagReco;
  std::unique_ptr<HistFactory> CutPlots;
  std::unique_ptr<MCPileupReweight> pileup_weights;
  JetId btag_medium, eta_cut, twojet, onejet, softjet, secondjet,wide_softjet;
  MuonId muid_cut, softMuon,isoMuon, hardMuon;
  ElectronId softElectron;
  TopJetId topjet, topjetid, heptopjetid,wjetId;
  std::unique_ptr<JetCleaner> jet_preclean;
  
  uhh2::Event::Handle<double> iso;
  uhh2::Event::Handle<double> trigger;
  
  std::unique_ptr<Selection> iso_sel;
  std::unique_ptr<OrSelection> trigger_sel;
};

OptSelModule::OptSelModule(Context& ctx){
  //get rid of jets that are outside the range of jet corrections
  jet_preclean.reset(new JetCleaner(ctx, PtEtaCut(5, 5)));

  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);

  wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 5.0));
  hardMuon = AndId<Muon>(MuonIDLoose(), PtEtaCut(55.0, 2.4));
  softMuon = AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(27.0, 2.4));
  isoMuon =  AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(27.0, 2.4),MuonIso(0.15));
  softElectron = AndId<Electron>(ElectronID_Spring16_veto_noIso, PtEtaCut(50.0, 2.5));
  softjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  topjet = PtEtaCut(150.0, 2.4); 
  topjetid = AndId<TopJet>(Type2TopTag(110,220, Type2TopTag::MassType::groomed,btag_medium),Tau32());
  onejet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(130.0, 2.4)); secondjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(40.0, 2.4)); 

  common.reset(new CommonModules());
  common->set_jet_id(wide_softjet);
  common->set_electron_id(softElectron);
  common->set_muon_id(softMuon);
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->set_HTjetid(softjet);
  common->init(ctx);


  iso =  ctx.declare_event_output<double>("IsoCriterion");
  trigger  =  ctx.declare_event_output<double>("trigger");
  iso_sel.reset(new NMuonSelection(1,1,isoMuon));
  trigger_sel.reset(new OrSelection());
  trigger_sel->add(make_unique<TriggerSelection>("HLT_Mu50_v*"));
  trigger_sel->add(make_unique<TriggerSelection>("HLT_MuTk50_v*"));

  OptTree.reset(new OptTreeModule(ctx));
 
  Reco.reset(new BprimeReco(ctx));
  Reco->set_jetRecoId(secondjet);
  TopTagReco.reset(new BprimeReco(ctx,"TopTagReco"));
  TopTagReco->set_topjetRecoId(topjetid);

  lepton.reset(new PrimaryLepton(ctx));
  chi2_combo.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::chi2_combo,"BprimeReco","Chi2Dis"));
  chi2_combo->set_emptyHyp(true);
  ttbar.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::ttbar,"BprimeReco","TTbarDis"));
  ttbar->set_emptyHyp(true);
  toptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::cmsTopTag,"TopTagReco","TopTagDis"));
  toptagDis->set_emptyHyp(true);

  CutPlots.reset(new HistFactory(ctx));
  CutPlots->setEffiHistName("CutsOpt");
  CutPlots->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*"),make_unique<TriggerSelection>("HLT_IsoMu24_v*"),make_unique<TriggerSelection>("HLT_IsoTkMu24_v*")),"muonTrigger");
  CutPlots->addSelection(make_unique<NElectronSelection>(0,0),"0_eleCut");
  CutPlots->addSelection(make_unique<NMuonSelection>(1,1,softMuon),"1_muonCut");
  CutPlots->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  CutPlots->addOrSelection(make_uvec(make_unique<TwoDCut>(0.4,40),make_unique<NMuonSelection>(1,1,isoMuon)),"Iso");
  CutPlots->addOrSelection(make_uvec(make_unique<NMuonSelection>(1,1,hardMuon),make_unique<NMuonSelection>(1,1,isoMuon)),"MuonCut");
  CutPlots->addSelection(make_unique<NJetSelection>(2,-1,secondjet),"50GeV_JetCut");
  CutPlots->addSelection(make_unique<NTopJetSelection>(1,-1,topjet),"150GeV_TopJetCut");
  //CutPlots->addSelection(make_unique<ForwardJetPtEtaCut>(1.5,40),"ForwardJetCut");
  CutPlots->addSelection(make_unique<METSelection>(40,ctx),"40GeV_METCut");
  CutPlots->addSelection(make_unique<HTLepSelection>(ctx,100),"100_HTLep");
  CutPlots->addHists("ElectronHists","ElectronHists");
  CutPlots->addHists("MuonHists","MuonHists");
  CutPlots->addHists("EventHists","EventHists");
  CutPlots->addHists("JetHists","JetHists");
  CutPlots->addHists("TopJetHists","TopJetHists");
  CutPlots->addHists("VLQGenHists","VLQGenHists");
  CutPlots->addHists("TopTag_JetHists",topjetid);
  CutPlots->addHists("BTagged_JetHists",btag_medium);
}

bool OptSelModule::process(Event & event){
  jet_preclean->process(event);
  if(!common->process(event)) return false;
  lepton->process(event);
  if(!CutPlots->passAndFill(event))return false;
  if(trigger_sel->passes(event))
    event.set(trigger,0.);
  else
    event.set(trigger,1.);
  if(iso_sel->passes(event))
    event.set(iso,1.);
  else
    event.set(iso,.0);
  Reco->massReco(event);
  ttbar->process(event);
  chi2_combo->process(event);
  TopTagReco->TopJetReco(event,2);
  toptagDis->process(event);
  OptTree->process(event);
  return true;
}
// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(OptSelModule)
