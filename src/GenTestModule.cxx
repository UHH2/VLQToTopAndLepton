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


class GenTestModule: public AnalysisModule {
public:

  explicit GenTestModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;
  std::unique_ptr<CommonModules> common, jetlepcleaning;
  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> bBprimeFactory;
  std::unique_ptr<HistFactory> muonFactory, softCleaningMuonFactory;
  std::unique_ptr<HistFactory> topWMuonFactory, wMuonFactory;
  std::unique_ptr<HistFactory> muonTrigger;
  std::unique_ptr<Selection> HiggsFilter, ZFilter;
  std::unique_ptr<Selection> genMttbar, test;
  AndSelection  channelSel;
  JetId btag_medium;
  ElectronId softElectron;
  MuonId muid_cut, softMuon;
  JetId secondjet, onejet, softjet, wide_softjet, ak4ForwardId, ak4CentralId;;
  TopJetId topjet,topjetid, hardtopjet;
  //std::unique_ptr<BprimeReco> Reco;
  std::unique_ptr<JetCleaner> jet_preclean;
  bool mttbar_sample=false;

  //jet uncertainties
  uhh2::Event::Handle<std::vector<Jet>> jet_jer_up, jet_jer_down, jet_jec_up, jet_jec_down;
  uhh2::Event::Handle<MET> met_jer_up, met_jer_down, met_jec_up, met_jec_down;
  std::vector<uhh2::Event::Handle<MET>>  met_handles;// = {met_jer_up, met_jer_down, met_jec_up, met_jec_down};
  std::string jercor_string="", jeccor_string=""; 
  std::unique_ptr<JetCorrector>  corrector_jec_up, corrector_jec_down;
  std::unique_ptr<JetResolutionSmearer> corrector_jer_up, corrector_jer_down;
  bool do_jer_unc = false;
  bool do_jec_unc = false;
  bool do_jet_uncer = false;
};

GenTestModule::GenTestModule(Context& ctx):channelSel(ctx){
  //set up the ttbar high mass samples
  mttbar_sample=false;
  jercor_string = ctx.get("jersmear_direction", "nominal");
  jeccor_string = ctx.get("jecsmear_direction", "nominal");
  if(jercor_string=="custom"){
    do_jer_unc = true;
    jet_jer_up   = ctx.declare_event_output<std::vector<Jet>>("jet_jer_up");
    jet_jer_down = ctx.declare_event_output<std::vector<Jet>>("jet_jer_down");
    met_jer_up   = ctx.declare_event_output<MET>("met_jer_up");
    met_jer_down = ctx.declare_event_output<MET>("met_jer_down");
    corrector_jer_up.reset(new JetResolutionSmearer(ctx,JERSmearing::SF_13TeV_2016,"jet_jer_up",1));
    corrector_jer_down.reset(new JetResolutionSmearer(ctx,JERSmearing::SF_13TeV_2016,"jet_jer_down",-1));
  }
  if(jeccor_string=="custom"){
    do_jec_unc = true;
    jet_jec_up   = ctx.declare_event_output<std::vector<Jet>>("jet_jec_up");
    jet_jec_down = ctx.declare_event_output<std::vector<Jet>>("jet_jec_down");
    met_jec_up   = ctx.declare_event_output<MET>("met_jec_up");
    met_jec_down = ctx.declare_event_output<MET>("met_jec_down");
    corrector_jec_up.reset(new JetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, {},"jet_jec_up","met_jec_up",1));
    corrector_jec_down.reset(new JetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, {},"jet_jec_down","met_jec_down",-1));
  }
  do_jet_uncer = do_jec_unc && do_jer_unc;
  if(do_jet_uncer)	
    met_handles = {met_jer_up, met_jer_down, met_jec_up, met_jec_down};



  string sample(ctx.get("dataset_version"));
  if(sample == "TTbar_Tune"){
    genMttbar.reset(new MttbarGenSelection(0,700));
    mttbar_sample=true;
  }

  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  lepton.reset(new PrimaryLepton(ctx));

  //put all variables here so that changes are applied to all cuts similar
  double delR_2D  = 0.4;
  double pTrel_2D = 25.;
  double MET_val = 50. ;
  double HTLep_val = 240.;
  double hardjetpt = 150.;
  double minjetpt = 30.;


  muid_cut = AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(55.0, 2.1)); //MuonIDTight()
  softMuon = AndId<Muon>(MuonIDLoose(), PtEtaCut(55.0, 2.1));
  softElectron = AndId<Electron>(ElectronID_Spring16_veto_noIso, PtEtaCut(115.0, 2.4));
  onejet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(130.0, 2.4)); 
  secondjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4)); 
  softjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(minjetpt, 3.0));
  wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(minjetpt, 5.0));
  topjet = PtEtaCut(150.0, 2.4); 
  hardtopjet = PtEtaCut(hardjetpt, 2.4); 
  topjetid = AndId<TopJet>(Type2TopTag(150,210, Type2TopTag::MassType::groomed,btag_medium),Tau32());
  ak4ForwardId = PtEtaCut(30.0,5,-1,2);
  ak4CentralId = PtEtaCut(30.0,2,-1,-1);

  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));

  //get rid of jets that are outside the range of jet corrections and do a first cleaning
  jet_preclean.reset(new JetCleaner(ctx, wide_softjet));
 
  common.reset(new CommonModules());
  if(!do_jet_uncer)common->set_jet_id(wide_softjet);
  common->set_electron_id(softElectron);
  common->set_muon_id(softMuon);
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->set_HTjetid(softjet);
  common->init(ctx);
 
  jetlepcleaning.reset(new CommonModules());
  //disable reweighting that was already done in common
  jetlepcleaning->disable_mclumiweight();
  jetlepcleaning->disable_mcpileupreweight();
  jetlepcleaning->disable_jec();
  jetlepcleaning->disable_jersmear();
  jetlepcleaning->disable_lumisel();
  jetlepcleaning->disable_metfilters();
  jetlepcleaning->disable_pvfilter();
  jetlepcleaning->disable_jetpfidfilter();
  //put all id such that it runs hopefully correctly
  jetlepcleaning->set_jet_id(wide_softjet);
  jetlepcleaning->set_electron_id(softElectron);
  jetlepcleaning->set_muon_id(softMuon);
  //jetlepcleaning->switch_jetlepcleaner();
  jetlepcleaning->switch_jetPtSorter();
  jetlepcleaning->set_HTjetid(softjet);
  jetlepcleaning->init(ctx);
  
  channelSel.add<NElectronSelection>("0Electrons",0,0);
  channelSel.add<NMuonSelection>("1Muon",1,1,muid_cut);

  muonFactory.reset(new HistFactory(ctx));
  muonFactory->setEffiHistName("muonEffis");
  muonFactory->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Mu50_v*"),make_unique<TriggerSelection>("HLT_TkMu50_v*")),"muonTrigger");
  //muonFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu50_v*"),"muonTrigger");
  muonFactory->addSelection(make_unique<NElectronSelection>(0,0,softElectron),"0_eleCut");
  muonFactory->addSelection(make_unique<NMuonSelection>(1,1,softMuon),"1_softMuonCut");
  muonFactory->addSelection(make_unique<NMuonSelection>(1,1,muid_cut),"1_muonCut");
  if(!do_jet_uncer)muonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  else muonFactory->addJetUncSelection(make_uvec(make_unique<NJetSelection>(1,-1,softjet,jet_jer_up),make_unique<NJetSelection>(1,-1,softjet,jet_jer_down),make_unique<NJetSelection>(1,-1,softjet,jet_jec_up),make_unique<NJetSelection>(1,-1,softjet,jet_jec_down)),met_handles,"30GeV_JetCut");
  if(!do_jet_uncer)muonFactory->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
  else muonFactory->addJetUncSelection(make_uvec(make_unique<TwoDCut>(ctx,"jet_jer_up",delR_2D,pTrel_2D),make_unique<TwoDCut>(ctx,"jet_jer_down",delR_2D,pTrel_2D),make_unique<TwoDCut>(ctx,"jet_jec_up",delR_2D,pTrel_2D),make_unique<TwoDCut>(ctx,"jet_jec_down",delR_2D,pTrel_2D)),met_handles,"2DCut");
  muonFactory->addAnalysisModule(move(jetlepcleaning),"JetLep_Cleaning");
  if(!do_jet_uncer)muonFactory->addSelection(make_unique<NJetSelection>(2,-1,secondjet),"50GeV_JetCut");
  else muonFactory->addJetUncSelection(make_uvec(make_unique<NJetSelection>(2,-1,secondjet,jet_jer_up),make_unique<NJetSelection>(2,-1,secondjet,jet_jer_down),make_unique<NJetSelection>(2,-1,secondjet,jet_jec_up),make_unique<NJetSelection>(2,-1,secondjet,jet_jec_down)),met_handles,"50GeV_JetCut");
  muonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),to_string((int)hardjetpt)+"GeV_TopJetCut");
  if(!do_jet_uncer)muonFactory->addSelection(make_unique<METSelection>(MET_val,ctx),to_string((int)MET_val)+"GeV_METCut");
  else muonFactory->addJetUncSelection(make_uvec(make_unique<METSelection>(MET_val,ctx,"met_jer_up"),make_unique<METSelection>(MET_val,ctx,"met_jer_down"),make_unique<METSelection>(MET_val,ctx,"met_jec_up"),make_unique<METSelection>(MET_val,ctx,"met_jec_down")),met_handles,to_string((int)MET_val)+"GeV_METCut");
  if(!do_jet_uncer)muonFactory->addSelection(make_unique<HTLepSelection>(ctx,HTLep_val),to_string((int)HTLep_val)+"GeV_HTLep");
  else muonFactory->addJetUncSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val,"met_jer_up"),make_unique<HTLepSelection>(ctx,HTLep_val,"met_jer_down"),make_unique<HTLepSelection>(ctx,HTLep_val,"met_jec_up"),make_unique<HTLepSelection>(ctx,HTLep_val,"met_jec_down")),met_handles,to_string((int)HTLep_val)+"GeV_HTLep");

  muonFactory->addHists("GenJetHists","muonChannel_GenJetHists");
  muonFactory->addHists("ElectronHists","muonChannel_ElectronHists");
  muonFactory->addHists("MuonHists","muonChannel_MuonHists");
  muonFactory->addHists("EventHists","muonChannel_EventHists");
  muonFactory->addHists("EventKinematicHists","muonChannel_EventKinematicsHists");
  muonFactory->addHists("JetHists","muonChannel_JetHists");
  muonFactory->addHists("TopJetHists","muonChannel_TopJetHists");
  muonFactory->addHists("VLQGenHists","muonChannel_VLQGenHists");
  muonFactory->addHists("LuminosityHists","muonChannel_LumiHists");
  muonFactory->addHists("Central_muonChannel_JetHists",ak4CentralId);
  muonFactory->addHists("Forward_muonChannel_JetHists",ak4ForwardId);
  muonFactory->addHists("50GeV_muonChannel_JetHists",secondjet);

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
  topWMuonFactory->addSelection(make_unique<NElectronSelection>(0,0,softElectron),"0_eleCut");
  topWMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,softMuon),"1_sofMuonCut");
  topWMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  topWMuonFactory->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,secondjet),"50GeV_JetCut");
  // topWMuonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),"300GeV_TopJetCut");
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
  wMuonFactory->addSelection(make_unique<NElectronSelection>(0,0,softElectron),"0_eleCut");
  wMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,softMuon),"1_softMuonCut");
  wMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  wMuonFactory->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,secondjet),"50GeV_JetCut");
  //wMuonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),"300GeV_TopJetCut");
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

  //muonTrigger.reset(new HistFactory(ctx,"triggerEffis.txt"));
  //muonTrigger.reset(new HistFactory(ctx));
  //muonTrigger->addSelection(make_unique<GenNSelection>(13,1,1,30,-1),"1_GenSel");
  //muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"muonTrigger");
  //muonTrigger->addHists("MuonHists","triggerChannel_MuonHists");

  cout<<"Mu Presel set up done"<<endl;
}

bool GenTestModule::process(Event & event){
  if(mttbar_sample)
    if(!genMttbar->passes(event))return false;
  jet_preclean->process(event);
  if(event.jets->size()==0) return false;
  if(do_jet_uncer){
    event.set(jet_jer_up, *event.jets);
    event.set(jet_jer_down, *event.jets);
    event.set(jet_jec_up, *event.jets);
    event.set(jet_jec_down, *event.jets);
    
    event.set(met_jec_up, *event.met);
    event.set(met_jec_down, *event.met);

    corrector_jer_up->process(event);
    corrector_jer_down->process(event);
    corrector_jec_up->process(event);
    corrector_jec_down->process(event);
  }
  if(!common->process(event)) return false;
  if(do_jet_uncer){
    event.set(met_jer_up, *event.met);
    event.set(met_jer_down, *event.met);
  }
 
  lepton->process(event);
  /*
  if(!event.isRealData){
    vlqGenHists->fill(event);  
    wMuonFactory->passAndFill(event);
    topWMuonFactory->passAndFill(event);
  }
  */
  bool muonChannel = muonFactory->passAndFill(event); 
  return muonChannel;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(GenTestModule)
