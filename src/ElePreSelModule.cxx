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

  AndSelection  channelSel;
  JetId btag_medium;
  ElectronId softElectron, eleId_cut;
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

ElePreSelModule::ElePreSelModule(Context& ctx):channelSel(ctx){
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
  //get rid of jets that are outside the range of jet corrections
  jet_preclean.reset(new JetCleaner(ctx, PtEtaCut(15, 5)));


  //put all variables here so that changes are applied to all cuts similar
  double delR_2D  = 0.4;
  double pTrel_2D = 25.;
  double MET_val = 60. ;
  double HTLep_val = 290.;
  double hardtopjetpt = 190;

  //Version  = ctx.get("dataset_version", "<not set>");
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  lepton.reset(new PrimaryLepton(ctx));

  muid_cut = AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(50.0, 2.1));//MuonIDTight()
  softMuon = AndId<Muon>(MuonIDLoose(), PtEtaCut(50.0, 2.1));
  softElectron = AndId<Electron>(ElectronID_Spring16_loose_noIso, PtEtaCut(110.0, 2.4));
  eleId_cut =  AndId<Electron>(ElectronID_Spring16_tight_noIso, PtEtaCut(118.0, 2.4));
  onejet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(250.0, 2.4)); secondjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4)); 
  softjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 5.0));
  topjet = PtEtaCut(150.0, 2.4); 
  hardtopjet = PtEtaCut(hardtopjetpt, 2.4); 
  topjetid = AndId<TopJet>(Type2TopTag(150,210, Type2TopTag::MassType::groomed,btag_medium),Tau32());
  ak4ForwardId = PtEtaCut(30.0,5,-1,2);
  ak4CentralId = PtEtaCut(30.0,2,-1,-1);

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
  jetlepcleaning->init(ctx);
  jetlepcleaning->set_HTjetid(softjet);



  channelSel.add<NElectronSelection>("0Electrons",1,1);
  channelSel.add<NMuonSelection>("1Muon",0,0,muid_cut);


  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists")); 

  eleFactory.reset(new HistFactory(ctx));
  eleFactory->setEffiHistName("eleEffis");
  eleFactory->addOrSelection(make_uvec(make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),make_unique<TriggerSelection>("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*")),"eleTrigger");
  eleFactory->addSelection(make_unique<NMuonSelection>(0,0,softMuon),"0_softMuonCut");
  eleFactory->addSelection(make_unique<NElectronSelection>(1,1,softElectron),"1_softeleCut");
  eleFactory->addSelection(make_unique<NElectronSelection>(1,1,eleId_cut),"1_eleCut");
  if(!do_jet_uncer)eleFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  else eleFactory->addJetUncSelection(make_uvec(make_unique<NJetSelection>(1,-1,softjet,jet_jer_up),make_unique<NJetSelection>(1,-1,softjet,jet_jer_down),make_unique<NJetSelection>(1,-1,softjet,jet_jec_up),make_unique<NJetSelection>(1,-1,softjet,jet_jec_down)),met_handles,"30GeV_JetCut");
  if(!do_jet_uncer)eleFactory->addSelection(make_unique<TwoDCut>(delR_2D,pTrel_2D),"2DCut");
  else eleFactory->addJetUncSelection(make_uvec(make_unique<TwoDCut>(ctx,"jet_jer_up",delR_2D,pTrel_2D),make_unique<TwoDCut>(ctx,"jet_jer_down",delR_2D,pTrel_2D),make_unique<TwoDCut>(ctx,"jet_jec_up",delR_2D,pTrel_2D),make_unique<TwoDCut>(ctx,"jet_jec_down",delR_2D,pTrel_2D)),met_handles,"2DCut");
  eleFactory->addAnalysisModule(move(jetlepcleaning),"JetLep_Cleaning");
  if(!do_jet_uncer)eleFactory->addSelection(make_unique<NJetSelection>(2,-1,secondjet),"50GeV_JetCut");
  else eleFactory->addJetUncSelection(make_uvec(make_unique<NJetSelection>(2,-1,secondjet,jet_jer_up),make_unique<NJetSelection>(2,-1,secondjet,jet_jer_down),make_unique<NJetSelection>(2,-1,secondjet,jet_jec_up),make_unique<NJetSelection>(2,-1,secondjet,jet_jec_down)),met_handles,"50GeV_JetCut");
   //eleFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"250GeV_JetCut");
  eleFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),"190GeV_TopJetCut");
  if(!do_jet_uncer)eleFactory->addSelection(make_unique<METSelection>(MET_val,ctx),to_string((int)MET_val)+"GeV_METCut");
  else eleFactory->addJetUncSelection(make_uvec(make_unique<METSelection>(MET_val,ctx,"met_jer_up"),make_unique<METSelection>(MET_val,ctx,"met_jer_down"),make_unique<METSelection>(MET_val,ctx,"met_jec_up"),make_unique<METSelection>(MET_val,ctx,"met_jec_down")),met_handles,to_string((int)MET_val)+"GeV_METCut");
  if(!do_jet_uncer)eleFactory->addSelection(make_unique<HTLepSelection>(ctx,HTLep_val),to_string((int)HTLep_val)+"GeV_HTLep");
  else eleFactory->addJetUncSelection(make_uvec(make_unique<HTLepSelection>(ctx,HTLep_val,"met_jer_up"),make_unique<HTLepSelection>(ctx,HTLep_val,"met_jer_down"),make_unique<HTLepSelection>(ctx,HTLep_val,"met_jec_up"),make_unique<HTLepSelection>(ctx,HTLep_val,"met_jec_down")),met_handles,to_string((int)HTLep_val)+"GeV_HTLep");
  //eleFactory->addSelection(make_unique<NJetSelection>(4,-1,softjet),"4_JetCut");

  //eleFactory->addAnalysisModule(make_unique<JetCleaner>(secondjet));
  //eleFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"130GeV_JetCut");
  //eleFactory->addSelection(make_unique<NTopJetSelection>(2,-1,topjet),"150GeV_2TopJetCut");
  //eleFactory->addSelection(make_unique<TopJetMassCut>(50),"50GeV_TopJetPruned");
  //eleFactory->addSelection(make_unique<ForwardJetPtEtaCut>(1.5,40),"ForwardJetCut");
  //eleFactory->addSelection(make_unique<NSubJetCut>(2),"2_SubJetCut");
  //eleFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<HTLepSelection>(ctx,200)),"200_HTLep_or_TopTag");
  //eleFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<NJetSelection>(4,-1,softjet)),"4_JetCut_or_TopTag");

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
  jet_preclean->process(event);
  if(event.jets->size()==0) return false;
  //cout<<"======================================="<<endl;
  //cout<<" module jet size "<<event.jets->size()<<endl;
  if(do_jet_uncer){
    event.set(jet_jer_up,*event.jets);
    event.set(jet_jer_down,*event.jets);
    event.set(jet_jec_up,*event.jets);
    event.set(jet_jec_down,*event.jets);
    
    event.set(met_jec_up,*event.met);
    event.set(met_jec_down,*event.met);

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
  if(!event.isRealData && 1==0){
    //event.weight *= 0.95;
    //bBprimeFactory->passAndFill(event);
    vlqGenHists->fill(event);  
    wMuonFactory->passAndFill(event);
    topWMuonFactory->passAndFill(event);
  }
  //Reco->massReco(event);
  //channelSel.passes(event);  
  //muonTrigger->passAndFill(event);
 
  bool eleChannel = eleFactory->passAndFill(event); 
  return eleChannel;
 
}

// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(ElePreSelModule)
