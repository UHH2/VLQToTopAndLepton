#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/TriggerSelection.h" 
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetCorrections.h" 
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"
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


class EleTriggerMeasurement: public AnalysisModule {
public:

  explicit EleTriggerMeasurement(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;
  std::unique_ptr<CommonModules> common, jetlepcleaning;
  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> lepton, ht;
  std::unique_ptr<HistFactory> muonTrigger, elesel;
  std::unique_ptr<Selection> genMttbar;
  
  ElectronId softElectron, eleId_cut,lowtriggerele;;
  MuonId muid_cut, softMuon, selMuon;
  JetId secondjet, onejet, hardtriggerjet, softjet, wide_softjet, ak4ForwardId, ak4CentralId, eta_cut;
  TopJetId topjet,topjetid, hardtopjet;
  std::unique_ptr<JetCleaner> jet_preclean;
  bool mttbar_sample=false;

  //jet uncertainties
  std::string jercor_string="", jeccor_string="";
  
  std::unique_ptr<JetCorrector>  corrector_jec_up, corrector_jec_down;
  std::unique_ptr<JetCleaner>  clean_jec_up, clean_jec_down, clean_jer_up, clean_jer_down;
  std::unique_ptr<JetResolutionSmearer> corrector_jer_up, corrector_jer_down;
  bool do_jer_unc = false;
  bool do_jec_unc = false;
  bool do_jet_uncer = false;
  bool run_muonid = false, run_muontrigger = false, run_eleid =false;

  std::unique_ptr<AnalysisModule> SF_muonID, SF_electronID, SF_muonTrigger;
  uhh2::Event::Handle<double> weight;
  uhh2::Event::Handle<bool> singleEletrigger, jetEletrigger, photontrigger, photontriggernoHE, ecalhttrigger;

  std::unique_ptr<TriggerSelection> singleEleTriggerSel, jetEleTriggerSel, photonTriggerSel, photonTriggerNoHESel, EcalHTTriggerSel;
  
};

EleTriggerMeasurement::EleTriggerMeasurement(Context& ctx){
  weight = ctx.declare_event_output<double>("weight");
  singleEletrigger = ctx.declare_event_output<bool>("singleEletrigger");
  jetEletrigger = ctx.declare_event_output<bool>("jetEletrigger");
  photontrigger = ctx.declare_event_output<bool>("photontrigger");
  photontriggernoHE = ctx.declare_event_output<bool>("photonnoHEtrigger");
  ecalhttrigger = ctx.declare_event_output<bool>("ecalhttrigger");  

  singleEleTriggerSel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
  jetEleTriggerSel.reset(new TriggerSelection("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*"));
  photonTriggerSel.reset(new TriggerSelection("HLT_Photon175_v*"));
  photonTriggerNoHESel.reset(new TriggerSelection("HLT_Photon300_NoHE_v*"));
  EcalHTTriggerSel.reset(new TriggerSelection("HLT_ECALHT800_v*"));
  
  //set up the ttbar high mass samples
  jercor_string = ctx.get("jersmear_direction", "nominal");
  jeccor_string = ctx.get("jecsmear_direction", "nominal");
  
  //Muon ScaleFactors
  if(!ctx.get("MounIDScaleFactors","").empty()){ 
    SF_muonID.reset(new MCMuonScaleFactor(ctx, ctx.get("MounIDScaleFactors"), "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tight", "nominal")); 
    run_muonid = true;
  }
  if(!ctx.get("MuonTriggerScaleFactors","").empty()){
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, ctx.get("MuonTriggerScaleFactors"), "IsoMu50_OR_IsoTkMu50_PtEtaBins", 1, "muonTrigger", "nominal")); 
    run_muontrigger = true;
  }
  //EleScaleFactors
  if(!ctx.get("EleScaleFactors","").empty()){
    SF_electronID.reset(new MCElecScaleFactor(ctx,ctx.get("EleScaleFactors"),1,"eleid","nominal" ));
    run_eleid =true;
  }
  string sample(ctx.get("dataset_version"));
  if(sample == "TTbar_Tune"){
    genMttbar.reset(new MttbarGenSelection(0,700));
    mttbar_sample=true;
  }                          

  //put all variables here so that changes are applied to all cuts similar
  double delR_2D  = 0.4;
  double pTrel_2D = 25.;
  double MET_val = 50. ;
  double HTLep_val = 240.;
  double hardjetpt = 150.;
  double minjetpt = 30.;

  lepton.reset(new PrimaryLepton(ctx));
  eta_cut = PtEtaCut(minjetpt,2.4);
  ht.reset(new HTCalc(ctx, eta_cut));
  
  wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(minjetpt, 5.0));
  eleId_cut =  AndId<Electron>(PtEtaCut(120.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);
  lowtriggerele =  AndId<Electron>(PtEtaCut(55.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);
  muid_cut = AndId<Muon>(MuonIDTight(), PtEtaCut(55.0, 2.4)); //MuonIDTight()
  selMuon = AndId<Muon>(MuonIDTight(), PtEtaCut(55.0, 2.4));
  softMuon = AndId<Muon>(PtEtaCut(40.0, 2.4),MuonIDTight());
  softElectron = AndId<Electron>(PtEtaCut(40.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);
  onejet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(130.0, 2.4)); 
  secondjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4)); 
  softjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(minjetpt, 3.0));
  topjet = PtEtaCut(150.0, 2.4); 
  hardtopjet = PtEtaCut(hardjetpt, 2.4); 
  ak4ForwardId = PtEtaCut(30.0,5,-1,2);
  ak4CentralId = PtEtaCut(30.0,2,-1,-1);
  hardtriggerjet=  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(185.0, 2.4));

  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));

  //get rid of jets that are outside the range of jet corrections and do a first cleaning
  jet_preclean.reset(new JetCleaner(ctx, wide_softjet));
   
  common.reset(new CommonModules());
  common->set_jet_id(wide_softjet);
  common->set_electron_id(softElectron);
  common->set_muon_id(selMuon);
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter(); 
  common->set_HTjetid(softjet);
  common->init(ctx);
  
  muonTrigger.reset(new HistFactory(ctx));
  muonTrigger->setEffiHistName("muonEffis");
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*")),"muonTrigger");
  muonTrigger->addSelection(make_unique<NMuonSelection>(1,-1,selMuon),"1_selMuonCut");
  muonTrigger->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  //muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")),"eleTrigger");
  muonTrigger->addSelection(make_unique<NElectronSelection>(1,1,softElectron),"1_eleCut");
  muonTrigger->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");

  muonTrigger->addHists("GenJetHists","trigger_presel_GenJetHists");
  muonTrigger->addHists("ElectronHists","trigger_presel_ElectronHists");
  muonTrigger->addHists("MuonHists","trigger_presel_MuonHists");
  muonTrigger->addHists("EventHists","trigger_presel_EventHists");
  muonTrigger->addHists("EventKinematicHists","trigger_presel_EventKinematicsHists");
  muonTrigger->addHists("JetHists","trigger_presel_JetHists");
  muonTrigger->addHists("TopJetHists","trigger_presel_TopJetHists");
  muonTrigger->addHists("VLQGenHists","trigger_presel_VLQGenHists");
  muonTrigger->addHists("LuminosityHists","trigger_presel_LumiHists");
  muonTrigger->addHists("Central_trigger_presel_JetHists",ak4CentralId);
  muonTrigger->addHists("Forward_trigger_presel_JetHists",ak4ForwardId);
  muonTrigger->addHists("50GeV_trigger_presel_JetHists",secondjet);
  muonTrigger->addHists("Lepton2DHist","trigger_presel_Lepton2D");

  std::unique_ptr<OrSelection> ele_combi_trigger;
  ele_combi_trigger.reset(new OrSelection());
  unique_ptr<AndSelection> singleEle(new AndSelection(ctx));
  singleEle->add("ele115",move(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")));
  singleEle->add("elept120",move(make_unique<NElectronSelection>(1,-1,eleId_cut)));  
  unique_ptr<AndSelection> jetEle(new AndSelection(ctx));
  jetEle->add("ele50_jet165pt",move(make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")));
  jetEle->add("ak4_185",move(make_unique<NJetSelection>(1,-1,hardtriggerjet)));
  jetEle->add("elept55",move(make_unique<NElectronSelection>(1,-1,lowtriggerele)));
  ele_combi_trigger->add(move(singleEle));
  ele_combi_trigger->add(move(jetEle));

  std::unique_ptr<OrSelection> ele_jet_combi;
  ele_jet_combi.reset(new OrSelection());
  unique_ptr<AndSelection> singleEle_notrig(new AndSelection(ctx));
  singleEle_notrig->add("elept120",move(make_unique<NElectronSelection>(1,-1,eleId_cut)));  
  unique_ptr<AndSelection> jetEle_notrig(new AndSelection(ctx));
  jetEle_notrig->add("ak4_185",move(make_unique<NJetSelection>(1,-1,hardtriggerjet)));
  jetEle_notrig->add("elept55",move(make_unique<NElectronSelection>(1,-1,lowtriggerele)));
  ele_jet_combi->add(move(singleEle_notrig));
  ele_jet_combi->add(move(jetEle_notrig));
  
  elesel.reset(new HistFactory(ctx));  
  elesel->addSelection(move(ele_combi_trigger),"eleTrigger_jet_electron_req");
  elesel->addSelection(move(ele_jet_combi),"jet_electron_trigreq");
  elesel->addSelection(make_unique<NJetSelection>(1,-1,hardtriggerjet),"jet185_cut");
  elesel->addSelection(make_unique<NElectronSelection>(1,1,lowtriggerele),"ele_cut");
  elesel->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")),"eleTrigger_combi");
  elesel->addAndSelection(make_uvec(make_unique<NJetSelection>(1,-1,hardtriggerjet),make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")),"eleTrigger_Ele50_PFJet165_jetcut");
  elesel->addAndSelection(make_uvec(make_unique<NElectronSelection>(1,1,lowtriggerele),make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")),"eleTrigger_Ele50_PFJet165_elecut");
  elesel->addAndSelection(make_uvec(make_unique<NElectronSelection>(1,1,eleId_cut),make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")),"eleTrigger_Ele115_elecut");
  elesel->addSelection(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),"eleTrigger_Ele115");
  elesel->addSelection(make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*"),"eleTrigger_Ele50_PFJet165");
  
 
  elesel->addHists("MuonHists","trigger_MuonHists");
  elesel->addHists("EventHists","trigger_EventHists");
  elesel->addHists("ElectronHists","trigger_ElectronHists");
  elesel->addHists("JetHists","trigger_JetHists");
  elesel->addHists("Lepton2DHist","trigger_Lepton2D");
}

bool EleTriggerMeasurement::process(Event & event){
  if(mttbar_sample)
    if(!genMttbar->passes(event))return false;
  jet_preclean->process(event);
  if(event.jets->size()==0) return false;
  if(!common->process(event)) return false;
  if(run_muontrigger) SF_muonTrigger->process(event);
  if(run_muonid) SF_muonID->process(event);
  if(run_eleid) SF_electronID->process(event); 
  lepton->process(event);
  ht->process(event);
  if(!muonTrigger->passAndFill(event)) return false;
  
  elesel->passAndFill(event,1);
  event.set(weight,event.weight);

  event.set(singleEletrigger,singleEleTriggerSel->passes(event));
  event.set(jetEletrigger,jetEleTriggerSel->passes(event));
  event.set(photontrigger,photonTriggerSel->passes(event));
  event.set(photontriggernoHE,photonTriggerNoHESel->passes(event));
  event.set(ecalhttrigger,EcalHTTriggerSel->passes(event));  

  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(EleTriggerMeasurement)
