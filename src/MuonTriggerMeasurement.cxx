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


class MuonTriggerMeasurement: public AnalysisModule {
public:

  explicit MuonTriggerMeasurement(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;
  std::unique_ptr<CommonModules> common, jetlepcleaning;
  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> muonTrigger, elesel;
  std::unique_ptr<Selection> genMttbar;
  
  ElectronId softElectron, eleId_cut, eleId_photon, lowtriggerele;;
  MuonId muid_loose, muid_tight, muid_highpt, softMuon, selMuon;
  JetId secondjet, onejet, hardtriggerjet, softjet, wide_softjet, ak4ForwardId, ak4CentralId;
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

  std::unique_ptr<AnalysisModule> SF_muonID, SF_electronID, SF_muonTrigger, SF_muonTrk, SF_ele;
  uhh2::Event::Handle<double> weight;
  uhh2::Event::Handle<bool> muontrigger, muonCombitrigger, highptmuon, tightmuon, loosemuon;

  std::unique_ptr<OrSelection> MuonTriggerSel , MuonCombiTriggerSel;
  std::unique_ptr<Selection> highptmuon_sel, tightmuon_sel, loosemuon_sel;
  
};

MuonTriggerMeasurement::MuonTriggerMeasurement(Context& ctx){
  weight = ctx.declare_event_output<double>("weight");
  muontrigger = ctx.declare_event_output<bool>("muontrigger");
  //mettrigger = ctx.declare_event_output<bool>("");

  muonCombitrigger = ctx.declare_event_output<bool>("muonCombitrigger");
  
  highptmuon = ctx.declare_event_output<bool>("highptmuon");
  tightmuon = ctx.declare_event_output<bool>("tightmuon");
  loosemuon = ctx.declare_event_output<bool>("loosemuon");
  
  //set up the ttbar high mass samples
  jercor_string = ctx.get("jersmear_direction", "nominal");
  jeccor_string = ctx.get("jecsmear_direction", "nominal");
  
  //Muon ScaleFactors
  if(!ctx.get("MounIDScaleFactors","").empty()){ 
    //SF_muonID.reset(new MCMuonScaleFactor(ctx, ctx.get("MounIDScaleFactors"), "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tight", "nominal")); 
    SF_muonID.reset(new MCMuonScaleFactor(ctx, ctx.get("MounIDScaleFactors"), "MC_NUM_HighPtID_DEN_genTracks_PAR_pt_eta", 1, "highpt",true, "nominal")); 
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

  lepton.reset(new PrimaryLepton(ctx));

  //put all variables here so that changes are applied to all cuts similar
  double delR_2D  = 0.4;
  double pTrel_2D = 25.;
  double MET_val = 50. ;
  double HTLep_val = 240.;
  double hardjetpt = 150.;
  double minjetpt = 30.;

  wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(minjetpt, 5.0));
  eleId_cut =  AndId<Electron>(PtEtaCut(120.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);
  eleId_photon=  AndId<Electron>(PtEtaCut(250.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);
  lowtriggerele =  AndId<Electron>(PtEtaCut(55.0, 2.4),ElectronID_MVAGeneralPurpose_Spring16_tight);
  muid_tight = AndId<Muon>(MuonIDTight(), PtEtaCut(20.0, 2.4));
  muid_loose = AndId<Muon>(MuonIDLoose(), PtEtaCut(20.0, 2.4));
  muid_highpt = AndId<Muon>(MuonIDHighPt(), PtEtaCut(20.0, 2.4));
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
  common->set_muon_id(PtEtaCut(20.0, 2.4));
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->set_HTjetid(softjet);
  common->init(ctx);

  highptmuon_sel.reset(new NMuonSelection(1,-1,muid_highpt));
  tightmuon_sel.reset(new NMuonSelection(1,-1,muid_tight));
  loosemuon_sel.reset(new NMuonSelection(1,-1,muid_loose));
  
  MuonTriggerSel.reset(new OrSelection());
  MuonTriggerSel->add(make_unique<TriggerSelection>("HLT_Mu50_v*"));
  MuonTriggerSel->add(make_unique<TriggerSelection>("HLT_TkMu50_v*"));
    
  MuonCombiTriggerSel.reset(new OrSelection());
  MuonCombiTriggerSel->add(make_unique<TriggerSelection>("HLT_Mu50_v*"));
  MuonCombiTriggerSel->add(make_unique<TriggerSelection>("HLT_TkMu50_v*"));
  MuonCombiTriggerSel->add(make_unique<TriggerSelection>("HLT_PFMETNoMu90_PFMHTNoMu90_v*"));

  muonTrigger.reset(new HistFactory(ctx));
  muonTrigger->setEffiHistName("muonEffis");
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*")),"muonTrigger");
  muonTrigger->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*"),make_unique<TriggerSelection>("HLT_PFMETNoMu90_PFMHTNoMu90_v*")),"MuonCombiTrigger");
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_PFMETNoMu90_PFMHTNoMu90_v*"),"PFMETNoMu90");
  muonTrigger->addSelection(make_unique<NMuonSelection>(1,-1,muid_highpt),"1_muonCut_highpt");
  muonTrigger->addSelection(make_unique<NMuonSelection>(1,-1,muid_tight),"1_muonCut_tight");
  muonTrigger->addHists("MuonHists","trigger_MuonHists");
  muonTrigger->addHists("EventHists","trigger_EventHists");
  muonTrigger->addHists("ElectronHists","trigger_ElectronHists");
  muonTrigger->addHists("JetHists","trigger_JetHists");
  muonTrigger->addHists("Lepton2DHist","trigger_Lepton2D");

  std::unique_ptr<OrSelection> ele_combi_trigger;
  ele_combi_trigger.reset(new OrSelection());
  unique_ptr<AndSelection> singleEle(new AndSelection(ctx));
  singleEle->add("ele115",move(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")));
  singleEle->add("elept120",move(make_unique<NElectronSelection>(1,1,eleId_cut)));  
  unique_ptr<AndSelection> jetEle(new AndSelection(ctx));
  jetEle->add("ele50_jet165pt",move(make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")));
  jetEle->add("ak4_185",move(make_unique<NJetSelection>(1,-1,hardtriggerjet)));
  jetEle->add("elept55",move(make_unique<NElectronSelection>(1,1,lowtriggerele)));
  unique_ptr<AndSelection> singlePhoton(new AndSelection(ctx));
  singlePhoton->add("photon175",move(make_unique<TriggerSelection>("HLT_Photon175_v*")));
  singlePhoton->add("elept250",move(make_unique<NElectronSelection>(1,1,eleId_photon)));  
  singlePhoton->add("ele50_jet165pt",move(make_unique<TriggerVeto>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")));
  singlePhoton->add("ele115",move(make_unique<TriggerVeto>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")));
  ele_combi_trigger->add(move(singlePhoton));
  ele_combi_trigger->add(move(singleEle));
  ele_combi_trigger->add(move(jetEle));
 
  elesel.reset(new HistFactory(ctx));
  elesel->addSelection(move(ele_combi_trigger),"eleTrigger_jet_electron_req");
  elesel->addSelection(make_unique<NElectronSelection>(1,-1,eleId_cut),"1_eleCut");
  elesel->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
  elesel->addSelection(make_unique<NMuonSelection>(1,-1),"1_muon");
  //elesel->addSelection(make_unique<NMuonSelection>(1,-1,muid_tight),"1_muonCut_tight");
  //elesel->addSelection(make_unique<NMuonSelection>(1,-1,muid_highpt),"1_muonCut_highpt");
  
  elesel->addHists("GenJetHists","trigger_presel_GenJetHists");
  elesel->addHists("ElectronHists","trigger_presel_ElectronHists");
  elesel->addHists("MuonHists","trigger_presel_MuonHists");
  elesel->addHists("EventHists","trigger_presel_EventHists");
  elesel->addHists("EventKinematicHists","trigger_presel_EventKinematicsHists");
  elesel->addHists("JetHists","trigger_presel_JetHists");
  elesel->addHists("TopJetHists","trigger_presel_TopJetHists");
  elesel->addHists("VLQGenHists","trigger_presel_VLQGenHists");
  elesel->addHists("LuminosityHists","trigger_presel_LumiHists");
  elesel->addHists("Central_trigger_presel_JetHists",ak4CentralId);
  elesel->addHists("Forward_trigger_presel_JetHists",ak4ForwardId);
  elesel->addHists("50GeV_trigger_presel_JetHists",secondjet);
  elesel->addHists("Lepton2DHist","trigger_presel_Lepton2D");
 
}

bool MuonTriggerMeasurement::process(Event & event){
  if(mttbar_sample)
    if(!genMttbar->passes(event))return false;
  jet_preclean->process(event);
  if(event.jets->size()==0) return false;
  if(!common->process(event)) return false;
  if(run_muontrigger) SF_muonTrigger->process(event);
  if(run_muonid) SF_muonID->process(event);
  if(run_eleid) SF_electronID->process(event); 
  lepton->process(event);
  
  if(!elesel->passAndFill(event)) return false;
  muonTrigger->passAndFill(event);
 
  event.set(weight,event.weight);
  event.set(muontrigger,MuonTriggerSel->passes(event));
  event.set(muonCombitrigger,MuonCombiTriggerSel->passes(event));
  event.set(highptmuon,highptmuon_sel->passes(event));
  event.set(tightmuon,tightmuon_sel->passes(event));
  event.set(loosemuon,loosemuon_sel->passes(event));

  return true;
  
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(MuonTriggerMeasurement)
