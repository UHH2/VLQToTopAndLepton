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



using namespace std;
using namespace uhh2;



class SelectionModule: public AnalysisModule {
public:

  explicit SelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  string Version;
  std::unique_ptr<BprimeRecoHists> RecoHists;
  std::unique_ptr<BprimeDiscriminator> chi2_combo, chi2_tlep, chi2_thad, ttbar, toptagDis, heptoptagDis;
  std::unique_ptr<AnalysisModule> ht, lepton;
  std::unique_ptr<BprimeReco> Reco, CMSTopTagReco, HEPTopTagReco;
  std::unique_ptr<BprimeRecoHists> Reco_wHad, Reco_wLep;
  std::unique_ptr<BprimeGen> Gen;
  std::unique_ptr<HistFactory> TagPlots;
  std::unique_ptr<HistFactory> Chi2Plots, CMSPlots, HEPPlots;
  //std::unique_ptr<HistFactory> TopTagBTag, TopTagNoBTag, NoTopTagBTag, NoTopTagNoBtag;
  std::unique_ptr<BprimeHypHists> chi2_Hists, ttbar_Hists, toptag_Hists, HEPtoptag_Hists;
  std::unique_ptr<Selection> topLep,topHad;
  std::unique_ptr<Selection> cmstoptagSel;
  std::unique_ptr<AndSelection> AntiCMSTopTag_AntiBTag, AntiCMSTopTag_BTag, AntiCMSTopTag_BTags2, AntiHEPTopTag_AntiBTag, AntiHEPTopTag_BTag, AntiHEPTopTag_BTags2;
  std::unique_ptr<JetHists> btag_jetHists;
  std::unique_ptr<Hists> jetHists_sortbyeta;
  std::unique_ptr<Hists> twoDjetHists_sortbyeta;
  std::unique_ptr<MCPileupReweight> pileup_weights;
  //std::unique_ptr<TopJetCleaner> topjetcleaner;
  JetId btag_medium;
  TopJetId topjetid;
  TopJetId heptopjetid;
  uhh2::Event::Handle<double> MCPileupWeight;
};

SelectionModule::SelectionModule(Context& ctx){
  pileup_weights.reset(new MCPileupReweight(ctx));
  //MCPileupWeight = ctx.declare_event_input<double>("MCPileupWeight", "MCPileupWeight"); //declare_in_out<float>("MCPileupWeight", "MCPileupWeight", ctx);
  topjetid = AndId<TopJet>(CMSTopTag(),Tau32());
  //heptopjetid = HEPTopTag();
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  //btag_medium = CSVBTag(CSVBTag::WP_TIGHT);
  Reco.reset(new BprimeReco(ctx));
  //Reco->set_jetRecoId(PtEtaCut(40.0,5));
  CMSTopTagReco.reset(new BprimeReco(ctx,"CMSTopTagReco"));
  CMSTopTagReco->set_topjetRecoId(topjetid);
  //HEPTopTagReco.reset(new BprimeReco(ctx,"HEPTopTagReco"));
  //HEPTopTagReco->set_topjetRecoId(heptopjetid);
  //HEPTopTagReco->set_topjetCollection(ctx,"patJetsHepTopTagCHSPacked");
  Gen.reset(new BprimeGen(ctx)); 
  ht.reset(new HTCalc(ctx));
  lepton.reset(new PrimaryLepton(ctx));

  jetHists_sortbyeta.reset(new JetHists(ctx,"jets_sortbyeta"));
  twoDjetHists_sortbyeta.reset(new VLQToTopAndLeptonHists(ctx,"twoDjets_sortbyeta"));
  btag_jetHists.reset(new JetHists(ctx,"btag_jetHists"));
  btag_jetHists->set_JetId(btag_medium);
  cmstoptagSel.reset(new NTopJetSelection(1,-1,topjetid));

  chi2_combo.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::chi2_combo));
  chi2_tlep.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::lepTop));
  chi2_thad.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::hadTop));
  ttbar.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::ttbar));
  toptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::cmsTopTag,"CMSTopTagReco"));
  //heptoptagDis.reset(new BprimeDiscriminator(ctx,BprimeDiscriminator::hepTopTag,"HEPTopTagReco"));
  ttbar_Hists.reset(new BprimeHypHists(ctx,"TTbarHists","DiscriminatorType_0"));
 
  TagPlots.reset(new HistFactory(ctx));
  TagPlots->setEffiHistName("Tags");
  TagPlots->addSelection(make_unique<NTopJetSelection>(1,-1,topjetid),"TopTag");
  TagPlots->addSelection(make_unique<NTopJetSelection>(0,0,topjetid),"AntiTopTag");
  //TagPlots->addSelection(make_unique<NTopJetSelection>(1,-1,heptopjetid),"HEPTopTag");
  //TagPlots->addSelection(make_unique<NTopJetSelection>(0,0,heptopjetid),"HEPAntiTopTag");
  TagPlots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"BTag");
  TagPlots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"2_BTags");
  TagPlots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"AntiBTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"AntiTopTag_AntiBTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"AntiTopTag_BTag");
  TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"AntiTopTag_2_BTag");
  //TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(0,0,btag_medium)),"AntiHEPTopTag_AntiBTag");
  //TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(1,1,btag_medium)),"AntiHEPTopTag_BTag");
  //TagPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"AntiHEPTopTag_2_BTag");
  TagPlots->addHists("ElectronHists","Tag_ElectronHists");
  TagPlots->addHists("MuonHists","Tag_MuonHists");
  TagPlots->addHists("EventHists","Tag_EventHists");
  TagPlots->addHists("JetHists","Tag_JetHists");
  TagPlots->addHists("TopJetHists","Tag_TopJetHists");
  TagPlots->addHists("VLQGenHists","Tag_VLQGenHists");

  Chi2Plots.reset(new HistFactory(ctx));
  Chi2Plots->setEffiHistName("Chi2Reco");
  Chi2Plots->addSelection(make_unique<PtRatioWTCut>(ctx,0.7,-1,"DiscriminatorType_1"),"Chi2_pTratioWTCut");
  //Chi2Plots->addSelection(make_unique<ChiSquareCut>(ctx,-1,45,"DiscriminatorType_0"),"Chi2_ttbarChi2Cut");
  Chi2Plots->addAndSelection(make_uvec(make_unique<ChiSquareCut>(ctx,20,0,"DiscriminatorType_1",12), make_unique<ChiSquareCut>(ctx,40,0,"DiscriminatorType_1",11)),"Chi2_BprimeChi2Cut");
  Chi2Plots->addSelection(make_unique<PTWhadCut>(ctx,200,-1,"DiscriminatorType_1"),"Chi2_PT_WhadCut");
  Chi2Plots->addSelection(make_unique<NTopJetSelection>(1,-1,topjetid),"Chi2_TopTag");
  Chi2Plots->addSelection(make_unique<NTopJetSelection>(0,0,topjetid),"Chi2_AntiTopTag");
  //Chi2Plots->addSelection(make_unique<NTopJetSelection>(1,-1,heptopjetid),"Chi2_HEPTopTag");
  //Chi2Plots->addSelection(make_unique<NTopJetSelection>(0,0,heptopjetid),"Chi2_AntiHEPTopTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"Chi2_BTag");
  Chi2Plots->addSelection(make_unique<NJetSelection>(2,2,btag_medium),"Chi2_2_BTags");
  Chi2Plots->addSelection(make_unique<NJetSelection>(3,3,btag_medium),"Chi2_3_BTags");
  Chi2Plots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"Chi2_2plus_BTags");
  Chi2Plots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"Chi2_AntiBTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"Chi2_AntiTopTag_AntiBTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"Chi2_AntiTopTag_BTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"Chi2_AntiTopTag_2_BTag");
  Chi2Plots->addAndSelection(make_uvec(make_unique<PTWhadCut>(ctx,200,-1,"DiscriminatorType_1"),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"Chi2_AntiTopTag_AntiBTag_WhadCut");
  Chi2Plots->addAndSelection(make_uvec(make_unique<PTWhadCut>(ctx,200,-1,"DiscriminatorType_1"),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"Chi2_AntiTopTag_BTag_WhadCut");
  Chi2Plots->addAndSelection(make_uvec(make_unique<PTWhadCut>(ctx,200,-1,"DiscriminatorType_1"),make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"Chi2_AntiTopTag_2_BTag_WhadCut");
  //Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(0,0,btag_medium)),"Chi2_AntiHEPTopTag_AntiBTag");
  //Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(1,1,btag_medium)),"Chi2_AntiHEPTopTag_BTag");
  //Chi2Plots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,heptopjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"Chi2_AntiHEPTopTag_2_BTag");
  Chi2Plots->addHists("ElectronHists","Chi2_ElectronHists");
  Chi2Plots->addHists("MuonHists","Chi2_MuonHists");
  Chi2Plots->addHists("EventHists","Chi2_EventHists");
  Chi2Plots->addHists("JetHists","Chi2_JetHists");
  Chi2Plots->addHists("Chi2_BJetHists",btag_medium);
  Chi2Plots->addHists("TopJetHists","Chi2_TopJetHists");
  Chi2Plots->addHists("Chi2_CMSTopTagJetHists",topjetid);
  Chi2Plots->addHists("VLQGenHists","Chi2_VLQGenHists");
  Chi2Plots->addHists("BprimeHypHists","Chi2_BprimeHypHists","DiscriminatorType_1");
  //Chi2Plots->addHists("BprimeHypHists","Chi2_tlep_BprimeHypHists","DiscriminatorType_2");
  //Chi2Plots->addHists("BprimeHypHists","Chi2_thad_BprimeHypHists","DiscriminatorType_3");

  CMSPlots.reset(new HistFactory(ctx));
  CMSPlots->setEffiHistName("CMSReco");
  CMSPlots->addSelection(make_unique<PtRatioWTCut>(ctx,0.7,-1,"DiscriminatorType_4"),"CMSReco_pTratioWTCut");
  CMSPlots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"CMSReco_BTag");
  CMSPlots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"CMSReco_2_BTags");
  CMSPlots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"CMSReco_AntiBTag");
  /*
  CMSPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"CMSReco_AntiTopTag_AntiBTag");
  CMSPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"CMSReco_AntiTopTag_BTag");
  CMSPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"CMSReco_AntiTopTag_2_BTag");*/
  CMSPlots->addHists("ElectronHists","CMSReco_ElectronHists");
  CMSPlots->addHists("MuonHists","CMSReco_MuonHists");
  CMSPlots->addHists("EventHists","CMSReco_EventHists");
  CMSPlots->addHists("JetHists","CMSReco_JetHists");
  CMSPlots->addHists("TopJetHists","CMSReco_TopJetHists");
  CMSPlots->addHists("CMSReco_CMSTopTagJetHists",topjetid);
  CMSPlots->addHists("VLQGenHists","CMSReco_VLQGenHists");
  CMSPlots->addHists("BprimeHypHists","CMSReco_BprimeHypHists","DiscriminatorType_4");

  /*
  HEPPlots.reset(new HistFactory(ctx));
  HEPPlots->setEffiHistName("HEPReco");
  HEPPlots->addSelection(make_unique<PtRatioWTCut>(ctx,0.7,-1,"DiscriminatorType_5"),"HEPReco_pTratioWTCut");
  HEPPlots->addSelection(make_unique<NJetSelection>(1,1,btag_medium),"HEPReco_BTag");
  HEPPlots->addSelection(make_unique<NJetSelection>(2,-1,btag_medium),"HEPReco_2_BTags");
  HEPPlots->addSelection(make_unique<NJetSelection>(0,0,btag_medium),"HEPReco_AntiBTag");
  HEPPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(0,0,btag_medium)),"HEPReco_AntiTopTag_AntiBTag");
  HEPPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(1,1,btag_medium)),"HEPReco_AntiTopTag_BTag");
  HEPPlots->addAndSelection(make_uvec(make_unique<NTopJetSelection>(0,0,topjetid),make_unique<NJetSelection>(2,-1,btag_medium)),"HEPReco_AntiTopTag_2_BTag");
  HEPPlots->addHists("ElectronHists","HEPReco_ElectronHists");
  HEPPlots->addHists("MuonHists","HEPReco_MuonHists");
  HEPPlots->addHists("EventHists","HEPReco_EventHists");
  HEPPlots->addHists("JetHists","HEPReco_JetHists");
  HEPPlots->addHists("TopJetHists","HEPReco_TopJetHists");
  HEPPlots->addHists("VLQGenHists","HEPReco_VLQGenHists");
  HEPPlots->addHists("BprimeHypHists","HEPReco_BprimeHypHists","DiscriminatorType_5");
  */
  vector<int> topLepIds {6,24,13};
  vector<int> topHadIds {6,24,-54321};
  topLep.reset(new GenFamilySelection(topLepIds,2));
  topHad.reset(new GenFamilySelection(topHadIds,2));;
  Reco_wHad.reset(new BprimeRecoHists(ctx, "Gen_wHad"));
  Reco_wLep.reset(new BprimeRecoHists(ctx, "Gen_wLep"));
}

bool SelectionModule::process(Event & event){
  pileup_weights->process(event);
  //if(!event.isRealData)event.weight *= event.get(MCPileupWeight);
  ht->process(event);
  lepton->process(event);
  if(!event.isRealData){
    Gen->process(event);  
    //lumi_weights->process(event);  
  }
  TagPlots->passAndFill(event,1);
  sort_by_eta(*event.jets);
  jetHists_sortbyeta->fill(event);
  btag_jetHists->fill(event);
  twoDjetHists_sortbyeta->fill(event);
  sort_by_pt(*event.jets);
  
  if(Reco->massReco(event)){
    if(ttbar->process(event))
      ttbar_Hists->fill(event);
    if(chi2_combo->process(event)){
      Chi2Plots->passAndFill(event,1);
    }
    if(!event.isRealData){
      if(topLep->passes(event))  Reco_wHad->fill(event);
      if(topHad->passes(event))  Reco_wLep->fill(event);
    }
  }
  /*
  if(HEPTopTagReco->TopJetReco(event,2)){
    if(heptoptagDis->process(event)){
      HEPPlots->passAndFill(event,1);
    }
  }
  */
  if(CMSTopTagReco->TopJetReco(event,2)){
    if(toptagDis->process(event)){
      CMSPlots->passAndFill(event,1);
    }
  }
  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SelectionModule)
