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


class TriggerSelModule: public AnalysisModule {
public:

  explicit TriggerSelModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> muonTrigger, eleTrigger, dataTrigger;
  std::unique_ptr<AndSelection> TopLepMuSel, TopLepEleLepSel, WLepMuSel, WLepEleSel, eleSel, muonSel;

  JetId btag_medium, jet, hardjet, jet180;
  MuonId muid;
  ElectronId eleid, lowtriggerele, eleid_lowpt;
  JetId twojet, onejet, lowjet;
   TopJetId topjet;
  bool first_event;

  //std::unique_ptr<BprimeReco> Reco;
};



TriggerSelModule::TriggerSelModule(Context& ctx){
  //Reco.reset(new BprimeReco(ctx)); 
  //Version  = ctx.get("dataset_version", "<not set>");
  first_event = false;
  jet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  common.reset(new CommonModules());
  common->disable_jec();
  common->disable_jersmear();
  common->disable_lumisel();
  common->disable_metfilters();
  common->disable_pvfilter();
  common->disable_jetpfidfilter();  
  common->set_jet_id(jet);
  common->set_electron_id(AndId<Electron>(ElectronID_MVAGeneralPurpose_Spring16_tight, PtEtaCut(120.0, 2.1)));
  common->set_muon_id(AndId<Muon>(MuonIDTight(),PtEtaCut(55.0, 2.1)));
  //common->switch_jetlepcleaner();
  common->set_HTjetid(jet);
  //common->switch_jetPtSorter();
  
  //common->disable_mclumiweight();
  //common->disable_mcpileupreweight();

  common->init(ctx);
  lepton.reset(new PrimaryLepton(ctx));
  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));

  //put all variables here so that changes are applied to all cuts similar
  double MET_val = 50. ;
  double HTLep_val = 150.;
  double ptrel = 25;
  double deltarmin =0.4;
  double st_cut = 2000;
  double met_cut= 250;
  jet180 = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(185.0, 2.4));
  hardjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(250.0, 2.4));
  
  muid  = AndId<Muon>(MuonIDTight(), PtEtaCut(55.0, 2.4));
  eleid = AndId<Electron>(PtEtaCut(120.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);
  eleid_lowpt = AndId<Electron>(PtEtaCut(20.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_loose);
 
  std::unique_ptr<OrSelection> muon_combi_trigger;
  muon_combi_trigger.reset(new OrSelection());
  unique_ptr<AndSelection> singleMuon(new AndSelection(ctx));
  singleMuon->add("HLT_Mu50",move(make_unique<TriggerSelection>("HLT_Mu50_v*")));
  unique_ptr<AndSelection> singleTkMuon(new AndSelection(ctx));
  singleTkMuon->add("HLT_TkMu50",move(make_unique<TriggerSelection>("HLT_TkMu50_v*")));
  unique_ptr<AndSelection> MetNoMu(new AndSelection(ctx));
  MetNoMu->add("HLT_METNoMu120",move(make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_v*")));
  MetNoMu->add("MET",move(make_unique<METSelection>(met_cut,ctx)));
  muon_combi_trigger->add(move(singleMuon));
  muon_combi_trigger->add(move(singleTkMuon));
  muon_combi_trigger->add(move(MetNoMu));

  muonTrigger.reset(new HistFactory(ctx));
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*")) ,"HLT_Mu50");
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*"),make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_v*")) ,"HLT_Mu50_OR_METNoMu");
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"),"HLT_METNoMu120");
  muonTrigger->addSelection(move(muon_combi_trigger),"HLT_Mu50_METNoMu120");
  muonTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"),make_unique<METSelection>(met_cut,ctx)),"HLT_METNoMu120_MET");
  
  muonTrigger->addHists("MuonHists","muon_trigger_MuonHists");
  muonTrigger->addHists("JetHists","muon_trigger_JetHists");
  muonTrigger->addHists("EventHists","muon_trigger_EventHists");
  muonTrigger->addHists("VLQGenHists","muon_trigger_VLQGenHists");

  lowtriggerele =  AndId<Electron>(PtEtaCut(55.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);

  std::vector<std::unique_ptr<OrSelection>> ele_combi_trigger_list;

  for(int i =0; i<8;i++){
    std::unique_ptr<OrSelection> ele_combi_trigger;
    ele_combi_trigger.reset(new OrSelection());
    unique_ptr<AndSelection> singleEle(new AndSelection(ctx));
    singleEle->add("ele115",move(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")));
    singleEle->add("elept120",move(make_unique<NElectronSelection>(1,-1,eleid)));
    if(i>4)
      singleEle->add("ecalht",move(make_unique<TriggerSelection>("HLT_ECALHT800_v*")));
    if(i>6)
      singleEle->add("ht",move(make_unique<STSelection>(ctx,st_cut)));
    unique_ptr<AndSelection> jetEle(new AndSelection(ctx));
    jetEle->add("ele50_jet165pt",move(make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")));
    jetEle->add("ak4_185",move(make_unique<NJetSelection>(1,-1,jet180)));
    jetEle->add("elept55",move(make_unique<NElectronSelection>(1,-1,lowtriggerele)));
    if(i>4)
      jetEle->add("ecalht",move(make_unique<TriggerSelection>("HLT_ECALHT800_v*")));
    //jetEle->add("ecalht",move(make_unique<TriggerSelection>("HLT_Photon175_v*")));
    if(i>6)
      jetEle->add("ht",move(make_unique<STSelection>(ctx,st_cut)));
    if(i>1 && i<5){
      unique_ptr<AndSelection> htEle(new AndSelection(ctx));
      htEle->add("ecalht",move(make_unique<TriggerSelection>("HLT_ECALHT800_v*")));
      htEle->add("ht",move(make_unique<STSelection>(ctx,st_cut)));
      ele_combi_trigger->add(move(htEle));
    }
    ele_combi_trigger->add(move(singleEle));
    ele_combi_trigger->add(move(jetEle));
    ele_combi_trigger_list.push_back(move(ele_combi_trigger));
  }

  eleTrigger.reset(new HistFactory(ctx));  
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_ECALHT800_v*"),"HLT_ECALHT800");
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"),"HLT_Ele105_CaloIdVT_GsfTrkIdT");	
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),"HLT_Ele115_CaloIdVT_GsfTrkIdT");	
  eleTrigger->addSelection(make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*"),"HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165");    
  eleTrigger->addAndSelection(make_uvec(move(ele_combi_trigger_list[0]),make_unique<NElectronSelection>(1,1,lowtriggerele)),"EleCombi_Trigger");
  eleTrigger->addAndSelection(make_uvec(move(ele_combi_trigger_list[2]),make_unique<NElectronSelection>(1,1,lowtriggerele)),"EleCombiHT800_Trigger");
  eleTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_ECALHT800_v*"),make_unique<STSelection>(ctx,st_cut)),"HT800_ST");
  //eleTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_ECALHT800_v*"),make_unique<STSelection>(ctx,1200)),"HT800_ST1200");
  eleTrigger->addHists("ElectronHists","ele_trigger_ElectronHists");
  eleTrigger->addHists("JetHists","ele_trigger_JetHists");
  eleTrigger->addHists("EventHists","ele_trigger_EventHists");
  eleTrigger->addHists("VLQGenHists","ele_trigger_VLQGenHists");


    
  std::unique_ptr<OrSelection> data_MuCombi;
  data_MuCombi.reset(new OrSelection());
  std::unique_ptr<AndSelection>data_Mu50trigger(new AndSelection(ctx));
  data_Mu50trigger->add("Mu50",make_unique<TriggerSelection>("HLT_Mu50_v*"));
  data_Mu50trigger->add("1muon",make_unique<NMuonSelection>(1,1,muid));
  std::unique_ptr<AndSelection>data_Mu50tkrtrigger(new AndSelection(ctx));
  data_Mu50tkrtrigger->add("TkMu50",make_unique<TriggerSelection>("HLT_TkMu50_v*"));
  data_Mu50tkrtrigger->add("1muon",make_unique<NMuonSelection>(1,1,muid));
  std::unique_ptr<AndSelection>data_METNoMutrigger(new AndSelection(ctx));
  data_METNoMutrigger->add("METNoMu",make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"));
  data_METNoMutrigger->add("MET",make_unique<METSelection>(met_cut));
  data_METNoMutrigger->add("1muon",make_unique<NMuonSelection>(1,1,muid));
  data_MuCombi->add(move(data_Mu50trigger));
  data_MuCombi->add(move(data_Mu50tkrtrigger));
  data_MuCombi->add(move(data_METNoMutrigger));
  
  dataTrigger.reset(new HistFactory(ctx));
  dataTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*")) ,"HLT_Mu50");
  dataTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*"),make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_v*")) ,"HLT_Mu50_OR_METNoMu120");
  dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_ECALHT800_v*"),"HLT_HT800");
  dataTrigger->addSelection(make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"),"HLT_PFMETNoMu120"); 
  dataTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"),make_unique<METSelection>(met_cut)),"PFMETNoMu120_MET");
  dataTrigger->addSelection(move(data_MuCombi),"MuCombi");
  
  dataTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_ECALHT800_v*"),make_unique<STSelection>(ctx,st_cut)),"HT800_ST");
  dataTrigger->addSelection(move(ele_combi_trigger_list[1]),"EleCombi_Sel");
  dataTrigger->addSelection(move(ele_combi_trigger_list[3]),"EleCombiHT800_Sel");
  dataTrigger->addSelection(move(ele_combi_trigger_list[5]),"HT800+EleCombi");
  dataTrigger->addSelection(move(ele_combi_trigger_list[7]),"HT800+ST+EleCombi");

  //dataTrigger->addAndSelection(make_uvec(make_unique<TriggerSelection>("HLT_ECALHT800_v*"),make_unique<STSelection>(ctx,st_cut),move(data_Mu50Comb)),"HT800+Mu50");
  //dataTrigger->addSelection(move(data_Mu50Comb),"HT800+Mu50");
	       
  dataTrigger->addHists("MuonHists","data_trigger_MuonHists");
  dataTrigger->addHists("ElectronHists","data_trigger_ElectronHists");
  dataTrigger->addHists("JetHists","data_trigger_JetHists");
  dataTrigger->addHists("EventHists","data_trigger_EventHists");

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


  eleSel.reset(new AndSelection(ctx)); muonSel.reset(new AndSelection(ctx));
  eleSel->add("electron",make_unique<NElectronSelection>(1,1)); muonSel->add("muon",make_unique<NMuonSelection>(1,1));
  eleSel->add("MET",make_unique<METSelection>(60));             muonSel->add("MET",make_unique<METSelection>(50));
  eleSel->add("HTLep",make_unique<HTLepSelection>(ctx,290));    muonSel->add("HTLep",make_unique<HTLepSelection>(ctx,250));
  eleSel->add("jet",make_unique<NJetSelection>(2,-1,jet));      muonSel->add("jet",make_unique<NJetSelection>(2,-1,jet));
  eleSel->add("2D",make_unique<TwoDCut>(deltarmin,ptrel));      muonSel->add("2D",make_unique<TwoDCut>(deltarmin,ptrel));
 
}

bool TriggerSelModule::process(Event & event){
  if(first_event){
     cout<<"Printing the trigger names for the first event"<<endl;
     for(auto item :  event.get_current_triggernames())
       cout<<item<<endl;
     first_event =false;
   }
  
  if(!common->process(event)) return false;
  lepton->process(event);
  if(!(eleSel->passes(event) || muonSel->passes(event)))
    return false;

  if(!event.isRealData){
    //if(topHadSel->passes(event) || topLepSel->passes(event))
    vlqGenHists->fill(event);
    if(TopLepMuSel->passes(event)||WLepMuSel->passes(event))
      muonTrigger->passAndFill(event,1);
    else if(TopLepEleLepSel->passes(event)||WLepEleSel->passes(event))
      eleTrigger->passAndFill(event,1);
    dataTrigger->passAndFill(event,1);
  }
  else{
    dataTrigger->passAndFill(event,1);
  }

  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TriggerSelModule)
