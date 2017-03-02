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
#include "UHH2/common/include/GenJetsHists.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TopPtReweight.h"

#include "UHH2/VLQToTopAndLepton/include/UncertaintyWeightsModule.h"

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
#include "UHH2/VLQToTopAndLepton/include/WTopJet.h"
//#include "UHH2/VLQToTopAndLepton/include/OptTreeModule.h"
#include "UHH2/VLQToTopAndLepton/include/JetReweight.h"
#include "UHH2/VLQToTopAndLepton/include/TopTagScalefactor.h"

#include "UHH2/VLQToTopAndLepton/include/WJetsReweight.h"


using namespace std;
using namespace uhh2;

class SelectionModule: public AnalysisModule {
public:

  explicit SelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  string Version;

  std::unique_ptr<CommonModules> common;
  std::unique_ptr<BprimeRecoHists> RecoHists;
  std::unique_ptr<BprimeDiscriminator> chi2_combo, chi2_tlep, chi2_thad, ttbar, cmstoptagDis, heptoptagDis, chi2_wtag, bestchi2_Dis;
  std::unique_ptr<AnalysisModule> ht, lepton, jetReweight;//, OptTree
  std::unique_ptr<BprimeReco> Reco, TopTagReco, HEPTopTagReco, WTagReco;
  std::unique_ptr<BprimeRecoHists> Reco_wHad, Reco_wLep;
  std::unique_ptr<BprimeGen> Gen;
  std::unique_ptr<HistFactory> TagPlots;
  std::unique_ptr<HistFactory> Chi2Plots, TopTagPlots, HEPPlots, WTagRecoPlots;
  //std::unique_ptr<HistFactory> TopTagBTag, TopTagNoBTag, NoTopTagBTag, NoTopTagNoBtag;
  std::unique_ptr<BprimeHypHists> chi2_Hists, ttbar_Hists, toptag_Hists, HEPtoptag_Hists;
  std::unique_ptr<Selection> topLep,topHad;
  std::unique_ptr<Selection> btagSel;
  std::unique_ptr<Selection> ttbar_chi2;
  std::unique_ptr<AndSelection> AntiTopTag_AntiBTag, AntiTopTag_BTag, AntiTopTag_BTags2, AntiHEPTopTag_AntiBTag, AntiHEPTopTag_BTag, AntiHEPTopTag_BTags2;
  std::unique_ptr<JetHists> btag_jetHists;
  std::unique_ptr<Hists> jetHists_sortbyeta, genjet_hists;
  std::unique_ptr<Hists> twoDjetHists_sortbyeta;
  //std::unique_ptr<TopJetCleaner> topjetcleaner;
  JetId subBtag, btag_medium,btag_tight,btag_loose, eta_cut, ak4ForwardId, ak4CentralId, jet;
  TopJetId topjetid,wjetId,heptopjetid;
  std::unique_ptr<Selection> hepselection;  
  std::unique_ptr<BTagMCEfficiencyHists> BTagEffiHists;
  std::unique_ptr<AnalysisModule> BTagScaleFactors,SF_muonID, SF_electronID, SF_muonTrigger;
  std::unique_ptr<AnalysisModule> pdf_scale_unc;
  std::unique_ptr<AnalysisModule> toptag_scale;
  std::unique_ptr<AnalysisModule> WJetsReweight_module;

  uhh2::Event::Handle<double> weight;
  uhh2::Event::Handle<double> numberofjets;

  string SF_muonID_variation ="nominal";
  //string SF_electronID_variation = "nominal";
  string BTag_variation ="central";
  string PU_variation ="central";

  bool run_muonid = false, run_muontrigger = false, run_eleid =false;

  std::unique_ptr<TopPtReweight> ttbar_reweight;

  uhh2::Event::Handle<FlavorParticle> primlep;

  int event_number = 0;
};

SelectionModule::SelectionModule(Context& ctx){
  weight = ctx.declare_event_output<double>("weight");
  numberofjets= ctx.declare_event_output<double>("numberofjets");
  primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  ttbar_reweight.reset(new TopPtReweight(ctx,0.159,-0.00141,"","weight_ttbar",true,0.9910819));
  pdf_scale_unc.reset(new UncertaintyWeightsModule(ctx));
  toptag_scale.reset(new TopTagScalefactor(ctx,"TopTagDis"));
  WJetsReweight_module.reset(new WJetsReweight(ctx));

  vector<int> topLepIds {6,24,13};
  vector<int> topHadIds {6,24,-54321};
  vector<int> wHadIds {24,-54321};
  vector<int> wLepIds {24,13};
  //topLep.reset(new GenFamilySelection(topLepIds,2));
  //topHad.reset(new GenFamilySelection(topHadIds,2));;
  //Reco_wHad.reset(new BprimeRecoHists(ctx, "Gen_wHad"));
  //Reco_wLep.reset(new BprimeRecoHists(ctx, "Gen_wLep"));

  //OptTree.reset(new OptTreeModule(ctx));
  jet = PtEtaCut(30,2.4);
  
  double forward_low = 2.0;
  double forward_upper = 5;
  double forwardJet_pt = 25;
  double drminForJet = 1000.0;
  double ForJetE = 250.;
  
  double central_low = 0.0;
  double central_upper = 2.0;
  double centralJet_pt = 25;
  double drminCenJet = 0.0;
  double CenJetE =0.0;

  //variations
  
  SF_muonID_variation = ctx.get("SF_muonID","nominal");
  BTag_variation =  ctx.get("BTag_variation","central");
  PU_variation = ctx.get("PU_variation","central");
  
  //common modules
  common.reset(new CommonModules());
  common->disable_jec();
  common->disable_jersmear();
  common->disable_jec();
  common->disable_lumisel();
  common->disable_metfilters();
  common->disable_pvfilter();
  common->disable_jetpfidfilter();
  common->init(ctx,PU_variation);

  //Muon ScaleFactors
  if(!ctx.get("MounIDScaleFactors","").empty()){ 
    SF_muonID.reset(new MCMuonScaleFactor(ctx, ctx.get("MounIDScaleFactors"), "MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta", 1, "mediummuon2016", "nominal")); 
    run_muonid = true;
  }
  if(!ctx.get("MuonTriggerScaleFactors","").empty()){
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, ctx.get("MuonTriggerScaleFactors"), "IsoMu50_OR_IsoTkMu50_PtEtaBins", 1, "muonTrigger", "nominal")); 
    run_muontrigger = true;
  }
  if(!ctx.get("EleScaleFactors","").empty()){ 
    SF_electronID.reset(new MCElecScaleFactor(ctx,ctx.get("EleScaleFactors"),1,"eleid","nominal" ));
    run_eleid =true;
  }

  //BTag Effi & Scale
  BTagEffiHists.reset(new BTagMCEfficiencyHists(ctx,"EffiHists/BTag",CSVBTag::WP_MEDIUM));
  BTagScaleFactors.reset(new MCBTagScaleFactor(ctx,CSVBTag::WP_MEDIUM,"jets",BTag_variation));

  wjetId = AndId<TopJet>(WMass(),Tau21(0.5));
  ak4ForwardId = PtEtaCut(30.0,5,-1,2);
  ak4CentralId = PtEtaCut(30.0,2,-1,-1);
  heptopjetid = HEPTopTagV2();
  hepselection.reset(new NTopJetSelection(1,-1,heptopjetid));
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  //subBtag = CSVBTag(0.79f);
  subBtag = CSVBTag(0.46f);
  btag_tight = CSVBTag(CSVBTag::WP_TIGHT);
  btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  //topjetid = AndId<TopJet>(Type2TopTag(110,210, Type2TopTag::MassType::groomed,subBtag),Tau32(0.54));
  topjetid = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed,subBtag),Tau32(0.5));
  Reco.reset(new BprimeReco(ctx));
  Reco->set_jetRecoId(jet);
  TopTagReco.reset(new BprimeReco(ctx,"TopTagReco"));
  TopTagReco->set_topjetRecoId(topjetid);

  Gen.reset(new BprimeGen(ctx)); 
  eta_cut = PtEtaCut(20.0,2.6);
  ht.reset(new HTCalc(ctx,eta_cut));
  lepton.reset(new PrimaryLepton(ctx));

  jetHists_sortbyeta.reset(new JetHists(ctx,"jets_sortbyeta"));
  genjet_hists.reset(new GenJetsHists(ctx,"genjets_hists"));

  twoDjetHists_sortbyeta.reset(new VLQToTopAndLeptonHists(ctx,"twoDjets_sortbyeta"));
  btag_jetHists.reset(new JetHists(ctx,"btag_jetHists"));
  btag_jetHists->set_JetId(btag_medium);
  btagSel.reset(new NJetSelection(1,-1,btag_medium));

  cmstoptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::cmsTopTag,"TopTagReco","TopTagDis"));
  cmstoptagDis->set_emptyHyp(true);
  ttbar.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::ttbar,"BprimeReco","TTbarDis"));
  ttbar->set_emptyHyp(true);
  chi2_combo.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::chi2_combo,"BprimeReco","Chi2Dis"));
  chi2_combo->set_emptyHyp(true);
  bestchi2_Dis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::chi_bestfit,"BprimeReco","BestFit"));
  bestchi2_Dis->set_emptyHyp(true);
  ttbar_Hists.reset(new BprimeHypHists(ctx,"TTbarHists","TTbarDis"));
  ttbar_chi2.reset(new ChiSquareCut(ctx,35,0,"TTbarDis"));

  TagPlots.reset(new HistFactory(ctx));
  TagPlots->setEffiHistName("Tags");
  TagPlots->addSelection(make_unique<NTopJetSelection>(1,-1,topjetid),"TopTag");
  TagPlots->addSelection(make_unique<NTopJetSelection>(0,0,topjetid),"AntiTopTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"AntiTopTag_AntiBTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"AntiTopTag_BTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"AntiTopTag_2_BTag");
  TagPlots->addSelection(make_unique<ForwardJetPtEtaCut>(ctx,forward_low,forward_upper,forwardJet_pt,-1,drminForJet,ForJetE),"ForwardJetEtaCut");
  TagPlots->addHists("ElectronHists","Tag_ElectronHists");
  TagPlots->addHists("MuonHists","Tag_MuonHists");
  TagPlots->addHists("EventHists","Tag_EventHists");
  TagPlots->addHists("JetHists","Tag_JetHists");
  TagPlots->addHists("TopJetHists","Tag_TopJetHists");
  TagPlots->addHists("VLQGenHists","Tag_VLQGenHists");
  TagPlots->addHists("EventKinematicHists","Tag_EventKinematicHists");
  TagPlots->addHists("TopTag_JetHists",topjetid);
  TagPlots->addHists("Forward_JetHists",ak4ForwardId);
  TagPlots->addHists("Central_JetHists",ak4CentralId);
  TagPlots->addHists("40GeV_JetHists",jet);
  TagPlots->addHists("BTaggedMedium_JetHists",btag_medium);
  TagPlots->addHists("BTaggedTight_JetHists",btag_tight);
  TagPlots->addHists("BTaggedLoose_JetHists",btag_loose);
  TagPlots->addHists("WTag_JetHists",wjetId);

  Chi2Plots.reset(new HistFactory(ctx));
  Chi2Plots->setEffiHistName("Chi2Reco");

  Chi2Plots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"Chi2_1_BTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"Chi2_2plus_BTags");
  Chi2Plots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"Chi2_AntiBTag");
  Chi2Plots->addSelection(make_unique<ForwardJetPtEtaCut>(ctx,0.,0.,0.,0.,-1,0.,"Chi2Dis"),"Chi2_NoForwardJet");


  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(1,1,btag_medium),make_unique<ForwardJetPtEtaCut> (ctx,forward_low,forward_upper,forwardJet_pt,-1,drminForJet,ForJetE,"Chi2Dis")),"Chi2_1_BTag_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,2,btag_medium),make_unique<ForwardJetPtEtaCut> (ctx,forward_low,forward_upper,forwardJet_pt,-1,drminForJet,ForJetE,"Chi2Dis")),"Chi2_2_BTags_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,-1,btag_medium),make_unique<ForwardJetPtEtaCut>(ctx,forward_low,forward_upper,forwardJet_pt,-1,drminForJet,ForJetE,"Chi2Dis")),"Chi2_2plus_BTags_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(0,0,btag_medium),make_unique<ForwardJetPtEtaCut> (ctx,forward_low,forward_upper,forwardJet_pt,-1,drminForJet,ForJetE,"Chi2Dis")),"Chi2_AntiBTag_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(1,1,btag_medium),make_unique<ForwardJetPtEtaCut> (ctx,central_low,central_upper,centralJet_pt,-1,drminCenJet,CenJetE,"Chi2Dis")),"Chi2_1_BTag_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,2,btag_medium),make_unique<ForwardJetPtEtaCut> (ctx,central_low,central_upper,centralJet_pt,-1,drminCenJet,CenJetE,"Chi2Dis")),"Chi2_2_BTags_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,-1,btag_medium),make_unique<ForwardJetPtEtaCut>(ctx,central_low,central_upper,centralJet_pt,-1,drminCenJet,CenJetE,"Chi2Dis")),"Chi2_2plus_BTags_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(0,0,btag_medium),make_unique<ForwardJetPtEtaCut> (ctx,central_low,central_upper,centralJet_pt,-1,drminCenJet,CenJetE,"Chi2Dis")),"Chi2_AntiBTag_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<GenFamilySelection>(topLepIds,2),make_unique<GenFamilySelection>(wHadIds,2)),"Chi2_TopLep_WHad");
  Chi2Plots->addAndSelection(make_uvec(make_unique<GenFamilySelection>(topHadIds,2),make_unique<GenFamilySelection>(wLepIds,2)),"Chi2_TopHad_WLep");

  Chi2Plots->addSelection(make_unique<ForwardJetPtEtaCut>(ctx,forward_low,forward_upper,forwardJet_pt,-1,drminForJet,ForJetE,"Chi2Dis"),"Chi2_Forward");
  Chi2Plots->addSelection(make_unique<ForwardJetPtEtaCut>(ctx,central_low,central_upper,centralJet_pt,-1,drminCenJet,CenJetE,"Chi2Dis"),"Chi2_Central");

  Chi2Plots->addHists("ElectronHists","Chi2_ElectronHists");
  Chi2Plots->addHists("MuonHists","Chi2_MuonHists");
  Chi2Plots->addHists("EventHists","Chi2_EventHists");
  Chi2Plots->addHists("JetHists","Chi2_JetHists");
  Chi2Plots->addHists("Chi2_BJetHists",btag_medium);
  Chi2Plots->addHists("Chi2_WTag",wjetId);
  Chi2Plots->addHists("TopJetHists","Chi2_TopJetHists");
  Chi2Plots->addHists("Chi2_TopTagJetHists",topjetid);
  Chi2Plots->addHists("VLQGenHists","Chi2_VLQGenHists");
  Chi2Plots->addHists("BprimeHypHists","Chi2_BprimeHypHists","Chi2Dis");
  //Chi2Plots->addHists("BprimeUncerHists","Chi2_BprimeUncerHists","Chi2Dis");
  Chi2Plots->addHists("EventKinematicHists","Chi2_EventKinematicHists_Dis","Chi2Dis");
  Chi2Plots->addHists("EventKinematicHists","Chi2_EventKinematicHists");

  TopTagPlots.reset(new HistFactory(ctx));
  //TopTagPlots->ScaleUncer();
  TopTagPlots->setEffiHistName("TopTagReco");
  TopTagPlots->addSelection(make_unique<ForwardJetPtEtaCut>(ctx,forward_low,forward_upper,forwardJet_pt,-1,drminForJet,ForJetE,"TopTagDis"),"TopTagReco_Forward");
  TopTagPlots->addSelection(make_unique<ForwardJetPtEtaCut>(ctx,central_low,central_upper,centralJet_pt,-1,drminCenJet,CenJetE,"TopTagDis"),"TopTagReco_Central");
  TopTagPlots->addAndSelection(make_uvec(make_unique<GenFamilySelection>(topLepIds,2),make_unique<GenFamilySelection>(wHadIds,2)),"TopTagReco_TopLep_WHad");
  TopTagPlots->addAndSelection(make_uvec(make_unique<GenFamilySelection>(topHadIds,2),make_unique<GenFamilySelection>(wLepIds,2)),"TopTagReco_TopHad_WLep");
  //TopTagPlots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"TopTagReco_2_BTag");
  TopTagPlots->addHists("ElectronHists","TopTagReco_ElectronHists");
  TopTagPlots->addHists("MuonHists","TopTagReco_MuonHists");
  TopTagPlots->addHists("EventHists","TopTagReco_EventHists");
  TopTagPlots->addHists("JetHists","TopTagReco_JetHists");
  TopTagPlots->addHists("TopJetHists","TopTagReco_TopJetHists");
  TopTagPlots->addHists("TopTagReco_TopTagJetHists",topjetid);
  TopTagPlots->addHists("VLQGenHists","TopTagReco_VLQGenHists");
  TopTagPlots->addHists("BprimeHypHists","TopTagReco_BprimeHypHists","TopTagDis");
  //TopTagPlots->addHists("BprimeUncerHists","TopTagReco_BprimeUncerHists","TopTagDis");
  TopTagPlots->addHists("EventKinematicHists","TopTagReco_EventKinematicHists_Dis","TopTagDis");
  TopTagPlots->addHists("EventKinematicHists","TopTagReco_EventKinematicHists");
}

bool SelectionModule::process(Event & event){ 
  WJetsReweight_module->process(event);
  common->process(event);
  ht->process(event);
  lepton->process(event);
  ttbar_reweight->process(event);
  Gen->process(event);  
  if(run_muontrigger) SF_muonTrigger->process(event);
  if(run_muonid) SF_muonID->process(event);
  if(run_eleid) SF_electronID->process(event);
  BTagEffiHists->fill(event);
  if(!event.isRealData){ 
    genjet_hists->fill(event);
  }
  TagPlots->passAndFill(event,1);
  sort_by_eta(*event.jets);
  jetHists_sortbyeta->fill(event);
  btag_jetHists->fill(event);
  twoDjetHists_sortbyeta->fill(event);
  sort_by_pt(*event.jets);

  bool reconstructed = false;
  //bool wreco = false; 

  bool toptagreco_bool = TopTagReco->TopJetReco(event,2);
  bool toptagdis_bool = cmstoptagDis->process(event);
  bool reco_bool = Reco->massReco(event);
  bool ttbardis_bool = ttbar->process(event);
  bool chi2dis_bool =  chi2_combo->process(event);
  bestchi2_Dis->process(event);

  if(toptagreco_bool && toptagdis_bool){
    reconstructed =true;
    TopTagPlots->passAndFill(event,1); 
  }

  BTagScaleFactors->process(event);
  if(reco_bool){
    if(ttbardis_bool)
      if(btagSel->passes(event) && ttbar_chi2->passes(event)) ttbar_Hists->fill(event); 
    if(chi2dis_bool && !reconstructed){
      Chi2Plots->passAndFill(event,1);
    }
  }
  toptag_scale->process(event);
  
  /*
   *
   * This is the version where I do not write out the tree and try to save some computing time
   * only considering certain possible reconstructions
   *
   */
  /*
  if(TopTagReco->TopJetReco(event,2)){
    if(cmstoptagDis->process(event)){
      reconstructed =true;
      TopTagPlots->passAndFill(event,1);
    }
  }
  if(Reco->massReco(event)){
    BTagScaleFactors->process(event);
    if(ttbar->process(event))
      if(btagSel->passes(event) && ttbar_chi2->passes(event)) ttbar_Hists->fill(event); 
    if(chi2_combo->process(event) && !reconstructed){
      Chi2Plots->passAndFill(event,1);
    }
  }
    */

  event.set(numberofjets,event.jets->size()); 
  event.set(weight,event.weight);
  pdf_scale_unc->process(event);
  

  return true;
}
// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SelectionModule)
