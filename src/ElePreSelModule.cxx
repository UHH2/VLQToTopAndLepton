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
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/TriggerSelection.h" 
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetCorrections.h" 
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/VLQToTopAndLepton/include/UncertaintyWeightsModule.h"

#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"
#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonHists.h"
#include "UHH2/VLQToTopAndLepton/include/GenSelection.h"
//#include "UHH2/VLQToTopAndLepton/include/HTSelection.h"
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"
#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"
#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"
#include "UHH2/VLQToTopAndLepton/include/Utils.h"
#include "UHH2/VLQToTopAndLepton/include/custom_jetcorr.h"

using namespace std;
using namespace uhh2;


class ElePreSelModule: public AnalysisModule {
public:

  explicit ElePreSelModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;
  std::unique_ptr<CommonModules> common, jetlepcleaning;
  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> bBprimeFactory;
  std::unique_ptr<HistFactory> eleFactory, softCleaningMuonFactory;
  std::unique_ptr<HistFactory> topWMuonFactory, wMuonFactory;
  std::unique_ptr<HistFactory> muonTrigger;
  std::unique_ptr<Selection> HiggsFilter, ZFilter;
  std::unique_ptr<Selection> genMttbar;

  JetId btag_medium;
  ElectronId softElectron, eleId_cut,lowtriggerele, highptele;
  MuonId muid_cut, softMuon;
  JetId secondjet, onejet, hardtriggerjet, softjet, wide_softjet, ak4ForwardId, ak4CentralId;
  TopJetId topjet,topjetid, hardtopjet;
  //std::unique_ptr<BprimeReco> Reco;
  std::unique_ptr<JetCleaner> jet_preclean;
  bool mttbar_sample=false;


  //jet uncertainties
  uhh2::Event::Handle<std::vector<Jet>> jet_jer_up, jet_jer_down, jet_jec_up, jet_jec_down;
  std::vector<uhh2::Event::Handle<std::vector<Jet>>> jet_collections;// = {jet_jer_up, jet_jer_down, jet_jec_up, jet_jec_down};
  uhh2::Event::Handle<MET> met_jer_up, met_jer_down, met_jec_up, met_jec_down;
  std::vector<uhh2::Event::Handle<MET>>  met_handles;// = {met_jer_up, met_jer_down, met_jec_up, met_jec_down};
  std::string jercor_string="", jeccor_string=""; 
  std::unique_ptr<JetCorrector>  corrector_jec_up, corrector_jec_down;
  std::unique_ptr<JetCleaner>  clean_jec_up, clean_jec_down, clean_jer_up, clean_jer_down;
  std::unique_ptr<JetResolutionSmearer> corrector_jer_up, corrector_jer_down;
  bool do_jer_unc = false;
  bool do_jec_unc = false;
  bool do_jet_uncer = false;
  
  string SF_muonID_variation ="nominal";
  //string SF_electronID_variation = "nominal";
  string BTag_variation ="central";
  string PU_variation ="central";
  std::unique_ptr<AnalysisModule> SF_muonID, SF_electronID, SF_muonTrigger;
  std::unique_ptr<custom_jetcorr> jet_uncer_provider;
  
  bool run_muonid = false, run_muontrigger = false, run_eleid =false;
};

ElePreSelModule::ElePreSelModule(Context& ctx){

  //put all variables here so that changes are applied to all cuts similar
  double delR_2D  = 0.4;
  double pTrel_2D = 25.;
  double MET_val = 60. ;
  double HTLep_val = 290.;
  double hardtopjetpt = 150;
  double minjetpt = 30.;
  string sample(ctx.get("dataset_version"));
  
  if(sample == "TTbar_Tune"){
    genMttbar.reset(new MttbarGenSelection(0,700));
    mttbar_sample=true;
  }  
  //get rid of jets that are outside the range of jet corrections
  jet_preclean.reset(new JetCleaner(ctx, PtEtaCut(15, 5)));
  if(jercor_string=="custom"){
    clean_jer_up.reset(new JetCleaner(ctx,wide_softjet,"jet_jer_up"));
    clean_jer_down.reset(new JetCleaner(ctx,wide_softjet,"jet_jer_down"));
  }
  if(jeccor_string=="custom"){
    clean_jec_up.reset(new JetCleaner(ctx,wide_softjet,"jet_jec_up"));
    clean_jec_down.reset(new JetCleaner(ctx,wide_softjet,"jet_jec_up"));
  }

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

  //Version  = ctx.get("dataset_version", "<not set>");
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  lepton.reset(new PrimaryLepton(ctx));

  muid_cut = AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(55.0, 2.1));//MuonIDTight()
  softMuon = AndId<Muon>(MuonIDLoose(), PtEtaCut(55.0, 2.1));
  softElectron = AndId<Electron>(ElectronID_Spring16_loose_noIso, PtEtaCut(55.0, 2.4));
  eleId_cut =  AndId<Electron>(ElectronID_Spring16_tight_noIso, PtEtaCut(120.0, 2.4));
  highptele =  AndId<Electron>(ElectronID_Spring16_tight_noIso, PtEtaCut(250.0, 2.4));
  lowtriggerele =  AndId<Electron>(ElectronID_Spring16_tight_noIso, PtEtaCut(55.0, 2.4));
  hardtriggerjet=  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(170.0, 2.4));
  onejet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(130.0, 2.4));
  secondjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4)); 
  softjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(minjetpt, 2.4));
  wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(minjetpt, 5.0));
  topjet = PtEtaCut(150.0, 2.4); 
  hardtopjet = PtEtaCut(hardtopjetpt, 2.4); 
  topjetid = AndId<TopJet>(Type2TopTag(150,210, Type2TopTag::MassType::groomed,btag_medium),Tau32());
  ak4ForwardId = PtEtaCut(30.0,5,-1,2);
  ak4CentralId = PtEtaCut(30.0,2,-1,-1);

  //get rid of jets that are outside the range of jet corrections and do a first cleaning
  jet_preclean.reset(new JetCleaner(ctx, wide_softjet));
  if(jercor_string=="custom"){
    clean_jer_up.reset(new JetCleaner(ctx,wide_softjet,"jet_jer_up"));
    clean_jer_down.reset(new JetCleaner(ctx,wide_softjet,"jet_jer_down"));
  }
  if(jeccor_string=="custom"){
    clean_jec_up.reset(new JetCleaner(ctx,wide_softjet,"jet_jec_up"));
    clean_jec_down.reset(new JetCleaner(ctx,wide_softjet,"jet_jec_up"));
  }

  common.reset(new CommonModules());
  if(!do_jet_uncer)common->set_jet_id(wide_softjet);
  common->set_electron_id(softElectron);
  common->set_muon_id(softMuon);
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->set_HTjetid(softjet);
  common->init(ctx);
  
  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));
  
  std::unique_ptr<OrSelection> trigger;
  trigger.reset(new OrSelection());
  unique_ptr<AndSelection> singleEle(new AndSelection(ctx));
  singleEle->add("ele115",move(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")));
  singleEle->add("elept120",move(make_unique<NElectronSelection>(1,1,eleId_cut)));
  unique_ptr<AndSelection> jetEle(new AndSelection(ctx));
  jetEle->add("ele50_jet165pt",move(make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")));
  jetEle->add("ak4_170",move(make_unique<NJetSelection>(1,-1,hardtriggerjet)));
  jetEle->add("elept55",move(make_unique<NElectronSelection>(1,1,lowtriggerele)));
  unique_ptr<AndSelection> photon(new AndSelection(ctx));
  photon->add("Photon175",move(make_unique<TriggerSelection>("HLT_Photon175_v*")));
  photon->add("ElePt250",move(make_unique<NElectronSelection>(1,1,highptele)));
  photon->add("ele50_jet165pt",move(make_unique<TriggerVeto>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")));
  photon->add("ele115",move(make_unique<TriggerVeto>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")));
  trigger->add(move(photon));
  //trigger->add(move(singleEle));
  //trigger->add(move(jetEle));
  
  eleFactory.reset(new HistFactory(ctx));
  eleFactory->setEffiHistName("eleEffis");
  eleFactory->addSelection(make_unique<NMuonSelection>(0,0,softMuon),"0_softMuonCut");
  eleFactory->addSelection(make_unique<NElectronSelection>(1,1,softElectron),"1_softeleCut");
  eleFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),to_string((int)hardtopjetpt)+"GeV_TopJetCut");

  if(!do_jet_uncer){
    eleFactory->addSelection(move(trigger),"eleTrigger");
    eleFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
    eleFactory->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
    eleFactory->addSelection(make_unique<NJetSelection>(2,-1,secondjet),"50GeV_JetCut");
    eleFactory->addSelection(make_unique<METSelection>(MET_val,ctx),to_string((int)MET_val)+"GeV_METCut");
    eleFactory->addSelection(make_unique<HTLepSelection>(ctx,HTLep_val),to_string((int)HTLep_val)+"GeV_HTLep");
  }


  if(do_jet_uncer){
    jet_uncer_provider.reset(new custom_jetcorr(ctx));
    met_handles = jet_uncer_provider->get_methandles();
    std::vector<std::unique_ptr<Selection>> widejet_sel;
    std::vector<std::unique_ptr<Selection>> twod_sel;	 
    std::vector<std::unique_ptr<AnalysisModule>> jetlep_sel;	 
    std::vector<std::unique_ptr<Selection>> hardjet_sel; 
    std::vector<std::unique_ptr<Selection>> met_sel;	 
    std::vector<std::unique_ptr<Selection>> htlep_sel;   

    std::vector<unique_ptr<Selection>> trigger_unc; 
    for(unsigned int i=0;i<4 && do_jet_uncer ;++i){
      unique_ptr<AndSelection> singleEle_unc(new AndSelection(ctx));
      singleEle_unc->add("ele115",move(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")));
      singleEle_unc->add("elept120",move(make_unique<NElectronSelection>(1,1,eleId_cut)));
      unique_ptr<AndSelection> jetEle_unc(new AndSelection(ctx));
      jetEle_unc->add("ele50_jet165pt",move(make_unique<TriggerSelection>("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*")));
      jetEle_unc->add("ak4_170",move(make_unique<NJetSelection>(1,-1,hardtriggerjet,jet_collections[i])));
      jetEle_unc->add("elept55",move(make_unique<NElectronSelection>(1,1,lowtriggerele)));
      unique_ptr<OrSelection> trigger_help;
      trigger_help.reset(new OrSelection);
      trigger_help->add(move(singleEle_unc));
      trigger_help->add(move(jetEle_unc));
      trigger_unc.push_back(move(trigger_help));
    }

    for(unsigned int i=0; i< met_handles.size();++i ){
      widejet_sel.push_back(make_unique<NJetSelection>(1,-1,softjet, jet_uncer_provider->get_jetHandle(i))); 
      twod_sel.push_back(make_unique<TwoDCut>(ctx,"jet_"+ jet_uncer_provider->get_namestring(i),delR_2D,pTrel_2D));	 
      jetlep_sel.push_back(make_unique<JetLeptonCleaner_by_KEYmatching>(ctx,JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC,"jet_"+ jet_uncer_provider->get_namestring(i)));	 
      hardjet_sel.push_back(make_unique<NJetSelection>(2,-1,secondjet,jet_uncer_provider->get_jetHandle(i)));     
      met_sel.push_back(make_unique<METSelection>(MET_val,ctx,"met_"+jet_uncer_provider->get_namestring(i)));	 
      htlep_sel.push_back(make_unique<HTLepSelection>(ctx,HTLep_val,"met_"+jet_uncer_provider->get_namestring(i)));   
    }
    eleFactory->addJetUncSelection(trigger_unc,met_handles,"trigger");
    eleFactory->addJetUncSelection(widejet_sel,met_handles,"30GeV_JetCut");
    eleFactory->addJetUncSelection(twod_sel,met_handles,"2DCut");
    eleFactory->addJetUncAnlysisModule(jetlep_sel,met_handles,"JetLep_Cleaning_jetcorr");
    eleFactory->addJetUncSelection(hardjet_sel,met_handles,"50GeV_JetCut");
    eleFactory->addJetUncSelection(met_sel,met_handles,to_string((int)MET_val)+"GeV_METCut");
    eleFactory->addJetUncSelection(htlep_sel,met_handles,to_string((int)HTLep_val)+"GeV_HTLep");
    return;
  }
  
  eleFactory->addHists("ElectronHists","eleChannel_ElectronHists");
  eleFactory->addHists("MuonHists","eleChannel_MuonHists");
  eleFactory->addHists("EventHists","eleChannel_EventHists");
  eleFactory->addHists("EventKinematicHists","eleChannel_EventKinematicsHists");
  eleFactory->addHists("JetHists","eleChannel_JetHists");
  eleFactory->addHists("TopJetHists","eleChannel_TopJetHists");
  eleFactory->addHists("VLQGenHists","eleChannel_VLQGenHists");
  eleFactory->addHists("LuminosityHists","eleChannel_LumiHists");
  eleFactory->addHists("Central_eleChannel_JetHists",ak4CentralId);
  eleFactory->addHists("Forward_eleChannel_JetHists",ak4ForwardId);
  eleFactory->addHists("50GeV_eleChannel_JetHists",secondjet);
  
  vector<int> bBprimeDecay {5, 6000007};

  bBprimeFactory.reset(new HistFactory(ctx));
  bBprimeFactory->addSelection(make_unique<GenFamilySelection>(bBprimeDecay,2),"bB'");


  vector<int> topLep {6,24,13};
  vector<int> topHad {6,24,-54321};
  vector<int> wHad {24,-54321};
  vector<int> wLep {24,13};
  topWMuonFactory.reset(new HistFactory(ctx));
  topWMuonFactory->setEffiHistName("topLep");
  topWMuonFactory->addSelection(make_unique<GenFamilySelection>(topLep,2),"topLep");
  topWMuonFactory->addSelection(make_unique<GenFamilySelection>(wHad,2),"wHad");
  topWMuonFactory->addSelection(make_unique<NMuonSelection>(0,0,softMuon),"0_sofMuonCut");
  topWMuonFactory->addSelection(make_unique<NElectronSelection>(1,1,softElectron),"1_softeleCut");
  topWMuonFactory->addSelection(make_unique<NElectronSelection>(1,1,eleId_cut),"1_eleCut");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  topWMuonFactory->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
  //topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"250GeV_JetCut");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,secondjet),"50GeV_JetCut");
  topWMuonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),"190GeV_TopJetCut");
  topWMuonFactory->addSelection(make_unique<METSelection>(MET_val,ctx),to_string((int)MET_val)+"GeV_METCut");
  topWMuonFactory->addSelection(make_unique<HTLepSelection>(ctx,HTLep_val),to_string((int)HTLep_val)+"GeV_HTLep");
  //topWMuonFactory->addSelection(make_unique<NJetSelection>(4,-1,softjet),"4_JetCut");
  

  topWMuonFactory->addHists("ElectronHists","topLep_ElectronHists");
  topWMuonFactory->addHists("MuonHists","topLep_MuonHists");
  topWMuonFactory->addHists("EventHists","topLep_EventHists");
  topWMuonFactory->addHists("EventKinematicHists","topLep_EventKinematicsHists");
  topWMuonFactory->addHists("JetHists","topLep_JetHists");
  topWMuonFactory->addHists("TopJetHists","topLep_TopJetHists");
  topWMuonFactory->addHists("VLQGenHists","topLep_VLQGenHists");
  topWMuonFactory->addHists("Central_topLep_JetHists",ak4CentralId);
  topWMuonFactory->addHists("Forward_topLep_JetHists",ak4ForwardId);
  topWMuonFactory->addHists("50GeV_topLep_JetHists",secondjet);

  wMuonFactory.reset(new HistFactory(ctx));
  wMuonFactory->setEffiHistName("topHad");	
  wMuonFactory->addSelection(make_unique<GenFamilySelection>(topHad,2),"topHad");
  wMuonFactory->addSelection(make_unique<GenFamilySelection>(wLep,2),"wLep");
  wMuonFactory->addSelection(make_unique<NMuonSelection>(0,0,softMuon),"0_softMuonCut");
  wMuonFactory->addSelection(make_unique<NElectronSelection>(1,1,softElectron),"1_softeleCut");
  wMuonFactory->addSelection(make_unique<NElectronSelection>(1,-1,eleId_cut),"1_eleCut");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  wMuonFactory->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
  //wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"250GeV_JetCut");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,secondjet),"50GeV_JetCut");
  wMuonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),"300GeV_TopJetCut");
  wMuonFactory->addSelection(make_unique<METSelection>(MET_val,ctx),to_string((int)MET_val)+"GeV_METCut");
  wMuonFactory->addSelection(make_unique<HTLepSelection>(ctx,HTLep_val),to_string((int)HTLep_val)+"GeV_HTLep");
  //wMuonFactory->addSelection(make_unique<NJetSelection>(4,-1,softjet),"4_JetCut");

  wMuonFactory->addHists("ElectronHists","topHad_ElectronHists");
  wMuonFactory->addHists("MuonHists","topHad_MuonHists");
  wMuonFactory->addHists("EventHists","topHad_EventHists");
  wMuonFactory->addHists("JetHists","topHad_JetHists");
  wMuonFactory->addHists("EventKinematicHists","topHad_EventKinematicsHists");
  wMuonFactory->addHists("TopJetHists","topHad_TopJetHists");
  wMuonFactory->addHists("VLQGenHists","topHad_VLQGenHists");
  wMuonFactory->addHists("Central_topHad_JetHists",ak4CentralId);
  wMuonFactory->addHists("Forward_topHad_JetHists",ak4ForwardId);
  wMuonFactory->addHists("50GeV_topHad_JetHists",secondjet);

  HiggsFilter.reset(new GenParticleFilter(25,0,0));
  ZFilter.reset(new GenParticleFilter(23,0,0));

}

bool ElePreSelModule::process(Event & event){
  if(mttbar_sample)
    if(!genMttbar->passes(event))return false;
  //if(genW_sample)
  //  if(!genW->passes(event))return false;
  jet_preclean->process(event);
  if(event.jets->size()==0) return false;
  if(do_jet_uncer){
    jet_uncer_provider->copy_jets(event);
    jet_uncer_provider->apply_corr(event);
    jet_uncer_provider->set_met(event);  
  } 
  if(!common->process(event)) return false;
  if(run_muontrigger) SF_muonTrigger->process(event);
  if(run_muonid) SF_muonID->process(event);
  if(run_eleid) SF_electronID->process(event);
  if(do_jet_uncer){
    jet_uncer_provider->clean(event);
  }
 
  lepton->process(event);
 
  bool electronChannel = eleFactory->passAndFill(event);
  return electronChannel;
}

// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(ElePreSelModule)
