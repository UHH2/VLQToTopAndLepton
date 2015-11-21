#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Selection.h"

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
#include "UHH2/VLQToTopAndLepton/include/WTopJet.h"
//#include "UHH2/VLQToTopAndLepton/include/OptTreeModule.h"

using namespace std;
using namespace uhh2;

class SelectionModule: public AnalysisModule {
public:

  explicit SelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  string Version;
  std::unique_ptr<BprimeRecoHists> RecoHists;
  std::unique_ptr<BprimeDiscriminator> chi2_combo, chi2_tlep, chi2_thad, ttbar, cmstoptagDis, heptoptagDis, chi2_wtag;
  std::unique_ptr<AnalysisModule> ht, lepton;//, OptTree
  std::unique_ptr<BprimeReco> Reco, CMSTopTagReco, HEPTopTagReco, WTagReco;
  std::unique_ptr<BprimeRecoHists> Reco_wHad, Reco_wLep;
  std::unique_ptr<BprimeGen> Gen;
  std::unique_ptr<HistFactory> TagPlots;
  std::unique_ptr<HistFactory> Chi2Plots, CMSPlots, HEPPlots, WTagRecoPlots;
  //std::unique_ptr<HistFactory> TopTagBTag, TopTagNoBTag, NoTopTagBTag, NoTopTagNoBtag;
  std::unique_ptr<BprimeHypHists> chi2_Hists, ttbar_Hists, toptag_Hists, HEPtoptag_Hists;
  std::unique_ptr<Selection> topLep,topHad;
  std::unique_ptr<Selection> btagSel;
  std::unique_ptr<Selection> ttbar_chi2;
  std::unique_ptr<AndSelection> AntiCMSTopTag_AntiBTag, AntiCMSTopTag_BTag, AntiCMSTopTag_BTags2, AntiHEPTopTag_AntiBTag, AntiHEPTopTag_BTag, AntiHEPTopTag_BTags2;
  std::unique_ptr<JetHists> btag_jetHists;
  std::unique_ptr<Hists> jetHists_sortbyeta;
  std::unique_ptr<Hists> twoDjetHists_sortbyeta;
  std::unique_ptr<MCPileupReweight> pileup_weights;
  //std::unique_ptr<TopJetCleaner> topjetcleaner;
  JetId btag_medium,btag_tight,btag_loose, eta_cut, ak4ForwardId, ak4CentralId, jet;
  TopJetId topjetid,wjetId,heptopjetid;
  std::unique_ptr<Selection> hepselection;  
};

SelectionModule::SelectionModule(Context& ctx){
  //OptTree.reset(new OptTreeModule(ctx));
  jet = PtEtaCut(40,2.4);
  pileup_weights.reset(new MCPileupReweight(ctx));
  topjetid = AndId<TopJet>(CMSTopTag(),Tau32());
  wjetId = AndId<TopJet>(WMass(),Tau21(0.5));
  ak4ForwardId = PtEtaCut(30.0,5,-1,2);
  ak4CentralId = PtEtaCut(30.0,2,-1,-1);
  heptopjetid = HEPTopTagV2();
  hepselection.reset(new NTopJetSelection(1,-1,heptopjetid));
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  btag_tight = CSVBTag(CSVBTag::WP_TIGHT);
  btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  Reco.reset(new BprimeReco(ctx));
  Reco->set_jetRecoId(jet);
  CMSTopTagReco.reset(new BprimeReco(ctx,"CMSTopTagReco"));
  CMSTopTagReco->set_topjetRecoId(topjetid);
  HEPTopTagReco.reset(new BprimeReco(ctx,"HEPTopTagReco"));
  HEPTopTagReco->set_topjetRecoId(heptopjetid);
  HEPTopTagReco->set_topjetCollection(ctx,"patJetsHepTopTagCHSPacked_daughters");
  Event::Handle<std::vector<TopJet>> heptopjets_handle = ctx.get_handle<std::vector<TopJet>>("patJetsHepTopTagCHSPacked_daughters");
  WTagReco.reset(new BprimeReco(ctx,"WTagReco"));
  WTagReco->set_wjetRecoId(wjetId);

  Gen.reset(new BprimeGen(ctx)); 
  eta_cut = PtEtaCut(20.0,2.6);
  ht.reset(new HTCalc(ctx,eta_cut));
  lepton.reset(new PrimaryLepton(ctx));

  jetHists_sortbyeta.reset(new JetHists(ctx,"jets_sortbyeta"));
  twoDjetHists_sortbyeta.reset(new VLQToTopAndLeptonHists(ctx,"twoDjets_sortbyeta"));
  btag_jetHists.reset(new JetHists(ctx,"btag_jetHists"));
  btag_jetHists->set_JetId(btag_medium);
  btagSel.reset(new NJetSelection(1,-1,btag_medium));

  chi2_combo.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::chi2_combo,"BprimeReco","Chi2Dis"));
  chi2_tlep.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::lepTop,"BprimeReco"));
  chi2_thad.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::hadTop,"BprimeReco"));
  ttbar.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::ttbar,"BprimeReco","TTbarDis"));
  chi2_wtag.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::wTag,"WTagReco","WTagDis"));
  cmstoptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::cmsTopTag,"CMSTopTagReco","CMSTopTagDis"));
  heptoptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::hepTopTag,"HEPTopTagReco","HEPTopTagDis"));
  ttbar_Hists.reset(new BprimeHypHists(ctx,"TTbarHists","TTbarDis"));
  ttbar_chi2.reset(new ChiSquareCut(ctx,35,0,"TTbarDis"));

  TagPlots.reset(new HistFactory(ctx));
  TagPlots->setEffiHistName("Tags");
  TagPlots->addSelection(make_unique<NTopJetSelection>(1,-1,topjetid),"TopTag");
  TagPlots->addSelection(make_unique<NTopJetSelection>(0,0,topjetid),"AntiTopTag");
  TagPlots->addSelection(make_unique<NTopJetSelection>(1,-1,heptopjetid,heptopjets_handle),"HEPTopTag");
  TagPlots->addSelection(make_unique<NTopJetSelection>(0,0,heptopjetid,heptopjets_handle),"HEPAntiTopTag");
  TagPlots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"BTag");
  TagPlots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"2_BTags");
  TagPlots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"AntiBTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"AntiTopTag_AntiBTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"AntiTopTag_BTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"AntiTopTag_2_BTag");
  TagPlots->addSelection(make_unique<ForwardJetPtEtaCut>(2.5,-1,30),"ForwardJetEta_2");
  //TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(0,0,btag_medium)),"AntiHEPTopTag_AntiBTag");
  //TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(1,1,btag_medium)),"AntiHEPTopTag_BTag");
  //TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"AntiHEPTopTag_2_BTag");
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

  Chi2Plots.reset(new HistFactory(ctx));
  Chi2Plots->setEffiHistName("Chi2Reco");
  //Chi2Plots->addSelection(make_unique<PtRatioWTCut>(ctx,0.7,-1,"DiscriminatorType_1"),"Chi2_pTratioWTCut");
  //Chi2Plots->addSelection(make_unique<ChiSquareCut>(ctx,-1,45,"DiscriminatorType_0"),"Chi2_ttbarChi2Cut");
  //Chi2Plots->addAndSelection(make_uvec(make_unique<ChiSquareCut>(ctx,20,0,"DiscriminatorType_1",12), make_unique<ChiSquareCut>(ctx,40,0,"DiscriminatorType_1",11)),"Chi2_BprimeChi2Cut");
  //Chi2Plots->addSelection(make_unique<PTWhadCut>(ctx,200,-1,"DiscriminatorType_1"),"Chi2_PT_WhadCut");
  Chi2Plots->addSelection(make_unique<NTopJetSelection>(1,-1,topjetid),"Chi2_TopTag");
  Chi2Plots->addSelection(make_unique<NTopJetSelection>(0,0,topjetid),"Chi2_AntiTopTag");
  //Chi2Plots->addSelection(make_unique<ForwardJetPtEtaCut>(2,-1,30),"ForwardJetEta_2");
  Chi2Plots->addSelection(make_unique<NTopJetSelection>(1,-1,heptopjetid,heptopjets_handle),"Chi2_HEPTopTag");
  Chi2Plots->addSelection(make_unique<NTopJetSelection>(0,0,heptopjetid,heptopjets_handle),"Chi2_AntiHEPTopTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(0,0,btag_loose),"Chi2_nolooseBTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(2,2,btag_tight),"Chi2_2_thighBTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(1,-1,btag_medium),"Chi2_BTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"Chi2_1_BTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(2,2,btag_medium),"Chi2_2_BTags");
  Chi2Plots->addSelection(make_unique<NJetSelection>(3,3,btag_medium),"Chi2_3_BTags");
  Chi2Plots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"Chi2_2plus_BTags");
  Chi2Plots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"Chi2_AntiBTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(1,1,btag_medium),make_unique<ForwardJetPtEtaCut>(2,5,30)),"Chi2_1_BTag_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,2,btag_medium),make_unique<ForwardJetPtEtaCut>(2,5,30)),"Chi2_2_BTags_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,-1,btag_medium),make_unique<ForwardJetPtEtaCut>(2,5,30)),"Chi2_2plus_BTags_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(0,0,btag_medium),make_unique<ForwardJetPtEtaCut>(2,5,30)),"Chi2_AntiBTag_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(1,1,btag_medium),make_unique<ForwardJetPtEtaCut>(0,2,30)),"Chi2_1_BTag_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,2,btag_medium),make_unique<ForwardJetPtEtaCut>(0,2,30)),"Chi2_2_BTags_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(2,-1,btag_medium),make_unique<ForwardJetPtEtaCut>(0,2,30)),"Chi2_2plus_BTags_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NJetSelection>(0,0,btag_medium),make_unique<ForwardJetPtEtaCut>(0,2,30)),"Chi2_AntiBTag_Central");
  Chi2Plots->addSelection(make_unique<ForwardJetPtEtaCut>(2,5,30),"Chi2_Forward");
  Chi2Plots->addSelection(make_unique<ForwardJetPtEtaCut>(0,2,30),"Chi2_Central");

  /*
  Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"Chi2_AntiTopTag_AntiBTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"Chi2_AntiTopTag_BTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"Chi2_AntiTopTag_2_BTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<ForwardJetPtEtaCut>(3,-1,30),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"Chi2_AntiTopTag_AntiBTag_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<ForwardJetPtEtaCut>(3,-1,30),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"Chi2_AntiTopTag_BTag_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<ForwardJetPtEtaCut>(3,-1,30),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"Chi2_AntiTopTag_2_BTag_Forward");
  Chi2Plots->addAndSelection(make_uvec(make_unique<ForwardJetPtEtaCut>(0,2,30),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"Chi2_AntiTopTag_AntiBTag_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<ForwardJetPtEtaCut>(0,2,30),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"Chi2_AntiTopTag_BTag_Central");
  Chi2Plots->addAndSelection(make_uvec(make_unique<ForwardJetPtEtaCut>(0,2,30),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"Chi2_AntiTopTag_2_BTag_Central");*/

  Chi2Plots->addHists("ElectronHists","Chi2_ElectronHists");
  Chi2Plots->addHists("MuonHists","Chi2_MuonHists");
  Chi2Plots->addHists("EventHists","Chi2_EventHists");
  Chi2Plots->addHists("JetHists","Chi2_JetHists");
  Chi2Plots->addHists("Chi2_BJetHists",btag_medium);
  Chi2Plots->addHists("TopJetHists","Chi2_TopJetHists");
  Chi2Plots->addHists("Chi2_CMSTopTagJetHists",topjetid);
  Chi2Plots->addHists("VLQGenHists","Chi2_VLQGenHists");
  Chi2Plots->addHists("BprimeHypHists","Chi2_BprimeHypHists","Chi2Dis");

  CMSPlots.reset(new HistFactory(ctx));
  CMSPlots->setEffiHistName("CMSReco");
  CMSPlots->addSelection(make_unique<ForwardJetPtEtaCut>(2,5,30),"CMSReco_Forward");
  CMSPlots->addSelection(make_unique<ForwardJetPtEtaCut>(0,2,30),"CMSReco_Central");
  CMSPlots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"CMSReco_2_BTag");
  CMSPlots->addHists("ElectronHists","CMSReco_ElectronHists");
  CMSPlots->addHists("MuonHists","CMSReco_MuonHists");
  CMSPlots->addHists("EventHists","CMSReco_EventHists");
  CMSPlots->addHists("JetHists","CMSReco_JetHists");
  CMSPlots->addHists("TopJetHists","CMSReco_TopJetHists");
  CMSPlots->addHists("CMSReco_CMSTopTagJetHists",topjetid);
  CMSPlots->addHists("VLQGenHists","CMSReco_VLQGenHists");
  CMSPlots->addHists("BprimeHypHists","CMSReco_BprimeHypHists","CMSTopTagDis");

  
  HEPPlots.reset(new HistFactory(ctx));
  HEPPlots->setEffiHistName("HEPReco");
  HEPPlots->addHists("ElectronHists","HEPReco_ElectronHists");
  HEPPlots->addHists("MuonHists","HEPReco_MuonHists");
  HEPPlots->addHists("EventHists","HEPReco_EventHists");
  HEPPlots->addHists("JetHists","HEPReco_JetHists");
  HEPPlots->addHists("TopJetHists","HEPReco_TopJetHists");
  HEPPlots->addHists("VLQGenHists","HEPReco_VLQGenHists");
  HEPPlots->addHists("BprimeHypHists","HEPReco_BprimeHypHists","HEPTopTagDis");

  WTagRecoPlots.reset(new HistFactory(ctx));
  WTagRecoPlots->setEffiHistName("WTagRecoReco");
  WTagRecoPlots->addSelection(make_unique<ForwardJetPtEtaCut>(2,5,30),"WTagReco_Forward");
  WTagRecoPlots->addSelection(make_unique<ForwardJetPtEtaCut>(0,2,30),"WTagReco_Central");
  WTagRecoPlots->addHists("ElectronHists","WTagReco_ElectronHists");
  WTagRecoPlots->addHists("MuonHists","WTagReco_MuonHists");
  WTagRecoPlots->addHists("EventHists","WTagReco_EventHists");
  WTagRecoPlots->addHists("JetHists","WTagReco_JetHists");
  WTagRecoPlots->addHists("TopJetHists","WTagReco_TopTagJetHists");
  WTagRecoPlots->addHists("WTagReco_WTagRecoJetHists",topjetid);
  WTagRecoPlots->addHists("VLQGenHists","WTagReco_VLQGenHists");
  WTagRecoPlots->addHists("BprimeHypHists","WTagReco_BprimeHypHists","WTagDis");


  vector<int> topLepIds {6,24,13};
  vector<int> topHadIds {6,24,-54321};
  topLep.reset(new GenFamilySelection(topLepIds,2));
  topHad.reset(new GenFamilySelection(topHadIds,2));;
  Reco_wHad.reset(new BprimeRecoHists(ctx, "Gen_wHad"));
  Reco_wLep.reset(new BprimeRecoHists(ctx, "Gen_wLep"));

}

bool SelectionModule::process(Event & event){
  pileup_weights->process(event);
  ht->process(event);
  lepton->process(event);
  if(!event.isRealData){
    Gen->process(event);  
  }
  TagPlots->passAndFill(event,1);
  sort_by_eta(*event.jets);
  jetHists_sortbyeta->fill(event);
  btag_jetHists->fill(event);
  twoDjetHists_sortbyeta->fill(event);
  sort_by_pt(*event.jets);

  if(HEPTopTagReco->TopJetReco(event,2)){
    if(heptoptagDis->process(event)){
      HEPPlots->passAndFill(event,1);
    }
  }
  bool reconstructed = false;
  //bool wreco = false; 
  if(CMSTopTagReco->TopJetReco(event,2)){
    if(cmstoptagDis->process(event)){
      reconstructed =true;
      CMSPlots->passAndFill(event,1);
    }
  }
  
  if(WTagReco->hadronicW(event,2) && !reconstructed){
    if(chi2_wtag->process(event))
      reconstructed =true;
      //wreco = true;
      WTagRecoPlots->passAndFill(event,1);
  }
  if(Reco->massReco(event)){
    if(ttbar->process(event))
      if(btagSel->passes(event) && ttbar_chi2->passes(event)) ttbar_Hists->fill(event); 
    if(chi2_combo->process(event) && !reconstructed){
      Chi2Plots->passAndFill(event,1);
    }
    if(!event.isRealData){
      if(topLep->passes(event))  Reco_wHad->fill(event);
      if(topHad->passes(event))  Reco_wLep->fill(event);
    }
  }
  return true;
}
// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SelectionModule)
