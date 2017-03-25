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
  std::unique_ptr<AndSelection> TopLepMuSel, TopLepEleLepSel, WLepMuSel, WLepEleSel;

  JetId btag_medium, jet, hardjet;
  MuonId muid_cut;
  JetId twojet, onejet, lowjet;
  TopJetId topjet;
  //std::unique_ptr<BprimeReco> Reco;
};



TriggerModule::TriggerModule(Context& ctx){
  //Reco.reset(new BprimeReco(ctx)); 
  //Version  = ctx.get("dataset_version", "<not set>");
  
  common.reset(new CommonModules());
  //common->set_jet_id(PtEtaCut(30.0,2.7));
  common->set_electron_id(AndId<Electron>(ElectronID_Spring16_medium_noIso, PtEtaCut(20.0, 2.1)));
  common->set_muon_id(AndId<Muon>(MuonIDMedium_ICHEP(),PtEtaCut(20.0, 2.1)));
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  //common->disable_mclumiweight();
  //common->disable_mcpileupreweight();
  common->init(ctx);
  lepton.reset(new PrimaryLepton(ctx));
  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));

  //put all variables here so that changes are applied to all cuts similar
  double MET_val = 50. ;
  double HTLep_val = 150.;
  double reliso_val = 0.15;
  double ptrel = 25;
  double deltarmin =0.4;
  jet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4));
  hardjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(250.0, 2.4));

  dataTrigger.reset(new HistFactory(ctx));
  dataTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*")) ,"HLT_Mu50");
  dataTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_IsoMu24_v*"),make_unique<TriggerSelection>("HLT_IsoTkMu24_v*")),"HLT_IsoMu24");
  dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*"),"HLT_Ele45_PFJet200_PFJet50");
  //dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele35_PFJet150_PFJet50_v*"),"HLT_Ele35_PFJet150_PFJet50");
  //dataTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_v*"),make_unique<TriggerSelection>("HLT_Ele45_PFJet200_PFJet50_v*")),"HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_AND_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50");
  dataTrigger->addHists("MuonHists","data_trigger_MuonHists");
  dataTrigger->addHists("ElectronHists","data_trigger_ElectronHists");
  dataTrigger->addHists("JetHists","data_trigger_JetHists");

  muonTrigger.reset(new HistFactory(ctx));
  //muonTrigger.reset(new HistFactory(ctx));
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*")) ,"HLT_Mu50");
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_IsoMu24_v*"),make_unique<TriggerSelection>("HLT_IsoTkMu24_v*")),"HLT_IsoMu24");			     		     
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v*"),"HLT_Mu40_eta2p1_PFJet200_PFJet50");
  muonTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TwoDCut>(deltarmin,ptrel)),"HLT_Mu50_2DIso");
  muonTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_IsoMu24_v*"),make_unique<RelIso>(reliso_val)),"HLT_IsoMu24_Iso");
  muonTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v*"),make_unique<NJetSelection>(2,-1,jet),make_unique<NJetSelection>(1,-1,hardjet)),"HLT_Mu40_eta2p1_PFJet200_PFJet50_req");
  muonTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v*"),make_unique<TwoDCut>(deltarmin,ptrel),make_unique<NJetSelection>(2,-1,jet),make_unique<NJetSelection>(1,-1,hardjet)),"HLT_Mu40_eta2p1_PFJet200_PFJet50_2DIso");
  muonTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val,ctx),make_unique<NJetSelection>(2,-1,jet)),"Sel");
  muonTrigger->addAndSelection(make_uvec(make_unique<RelIso>(reliso_val),make_unique<TwoDCut>(deltarmin,ptrel),make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val,ctx),make_unique<NJetSelection>(2,-1,jet)),"Sel_2D_Iso");
  muonTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val,ctx),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_Mu50_v*")),"HLT_Mu50_Sel");
  muonTrigger->addAndSelection(make_uvec(make_unique<TwoDCut>(deltarmin,ptrel),make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val,ctx),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_Mu50_v*")),"HLT_Mu50_Sel_2DCut");
  muonTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_TkMu50_v*")),"HLT_TkMu50_Sel");
  muonTrigger->addAndSelection(make_uvec(make_unique<RelIso>(reliso_val),make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_IsoMu24_v*")),"HLT_IsoMu24_Sel_Iso");  		       
  muonTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_IsoMu24_v*")),"HLT_IsoMu24_Sel");  		     

  muonTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v*")),"HLT_Mu40_eta2p1_PFJet200_PFJet50_Sel");
  muonTrigger->addHists("MuonHists","muon_trigger_MuonHists");
  muonTrigger->addHists("JetHists","muon_trigger_JetHists");
  muonTrigger->addHists("VLQGenHists","muon_trigger_VLQGenHists");

  eleTrigger.reset(new HistFactory(ctx));
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele27_eta2p1_WPTight_Gsf_v*"),"HLT_Ele27_eta2p1");
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*"),"HLT_Ele45_PFJet200_PFJet50");	
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"),"HLT_Ele105_CaloIdVT_GsfTrkIdT");	
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),"HLT_Ele115_CaloIdVT_GsfTrkIdT");	

  eleTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele27_eta2p1_WPTight_Gsf_v*"),make_unique<RelIso>(reliso_val)),"HLT_Ele27_eta2p1_Iso");
  eleTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*"),make_unique<NJetSelection>(2,-1,jet),make_unique<NJetSelection>(1,-1,hardjet)),"HLT_Ele45_PFJet200_PFJet50_req");	
  eleTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"),make_unique<TwoDCut>(deltarmin,ptrel)),"HLT_Ele105_CaloIdVT_GsfTrkIdT_2DIso");	
  eleTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),make_unique<TwoDCut>(deltarmin,ptrel)),"HLT_Ele115_CaloIdVT_GsfTrkIdT_2DIso");	
  eleTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet)),"Sel");	

  eleTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_Ele27_eta2p1_WPTight_Gsf_v*")),"HLT_Ele27_eta2p1_Sel");
  eleTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*")),"HLT_Ele45_PFJet200_PFJet50_Sel");	
  eleTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*")),"HLT_Ele105_CaloIdVT_GsfTrkIdT_Sel");
  eleTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(2,-1,jet),make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")),"HLT_Ele115_CaloIdVT_GsfTrkIdT_Sel");
  eleTrigger->addAndSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val),make_unique<METSelection>(MET_val),make_unique<NJetSelection>(1,-1,hardjet),make_unique<TriggerSelection>("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*")),"HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_Sel_with250GeVJet");
  eleTrigger->addHists("ElectronHists","ele_trigger_ElectronHists");
  eleTrigger->addHists("JetHists","ele_trigger_JetHists");
  eleTrigger->addHists("VLQGenHists","ele_trigger_VLQGenHists");
  
  vector<int> topLepMu {6,24,13};
  vector<int> topLepEle {6,24,11};
  vector<int> topHad {6,24,-54321};
  vector<int> wHad {24,-54321};
  vector<int> wLepMu {24,13};
  vector<int> wLepEle {24,11};

  TopLepMuSel.reset(new AndSelection(ctx));
  TopLepEleLepSel.reset(new AndSelection(ctx));
  WLepMuSel.reset(new AndSelection(ctx));
  WLepEleSel.reset(new AndSelection(ctx));

  TopLepMuSel->add("TopLepMuSel_lep",make_unique<GenFamilySelection>(topLepMu,2));
  TopLepMuSel->add("TopLepMuSel_had",make_unique<GenFamilySelection>(wHad,2));

  TopLepEleLepSel->add("TopLepEleLepSel_lep",make_unique<GenFamilySelection>(topLepEle,2));
  TopLepEleLepSel->add("TopLepEleLepSel_had",make_unique<GenFamilySelection>(wHad,2));

  WLepMuSel->add("WLepMuSel_lep",make_unique<GenFamilySelection>(wLepMu,2));
  WLepMuSel->add("WLepMuSel_had",make_unique<GenFamilySelection>(topHad,2));

  WLepEleSel->add("WLepEleSel_lep",make_unique<GenFamilySelection>(wLepEle,2));
  WLepEleSel->add("WLepEleSel_had",make_unique<GenFamilySelection>(topHad,2));

}

bool TriggerModule::process(Event & event){
  if(!common->process(event)) return false;
  lepton->process(event);
  if(!event.isRealData){
    //if(topHadSel->passes(event) || topLepSel->passes(event))
    vlqGenHists->fill(event);
    if(TopLepMuSel->passes(event)||WLepMuSel->passes(event))
      muonTrigger->passAndFill(event,1);
    else if(TopLepEleLepSel->passes(event)||WLepEleSel->passes(event))
      eleTrigger->passAndFill(event,1);
  }
  else{
    dataTrigger->passAndFill(event,1);
  }

  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TriggerModule)
