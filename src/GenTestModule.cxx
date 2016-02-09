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


class GenTestModule: public AnalysisModule {
public:

  explicit GenTestModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> bBprimeFactory;
  std::unique_ptr<HistFactory> muonFactory, softCleaningMuonFactory;
  std::unique_ptr<HistFactory> topWMuonFactory, wMuonFactory;
  std::unique_ptr<HistFactory> muonTrigger;
  std::unique_ptr<Selection> HiggsFilter, ZFilter;
  AndSelection  channelSel;
  JetId btag_medium;
  ElectronId softElectron;
  MuonId muid_cut, softMuon;
  JetId secondjet, onejet, softjet, wide_softjet, ak4ForwardId, ak4CentralId;;
  TopJetId topjet,topjetid, hardtopjet;
  //std::unique_ptr<BprimeReco> Reco;
};

GenTestModule::GenTestModule(Context& ctx):channelSel(ctx){
  //Version  = ctx.get("dataset_version", "<not set>");
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  lepton.reset(new PrimaryLepton(ctx));

  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));
  //muonTrigger.reset(new HistFactory(ctx,"triggerEffis.txt"));
  //muonTrigger.reset(new HistFactory(ctx));
  //muonTrigger->addSelection(make_unique<GenNSelection>(13,1,1,30,-1),"1_GenSel");
  //muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"muonTrigger");
  //muonTrigger->addHists("MuonHists","triggerChannel_MuonHists");
  muid_cut = AndId<Muon>(MuonIDTight(), PtEtaCut(50.0, 2.1));
  softMuon = AndId<Muon>(MuonIDLoose(), PtEtaCut(50.0, 2.1));
  softElectron = AndId<Electron>(ElectronID_Spring15_25ns_loose_noIso, PtEtaCut(50.0, 2.5));
  onejet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(130.0, 2.4)); secondjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(50.0, 2.4)); 
  softjet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 5.0));
  topjet = PtEtaCut(150.0, 2.4); 
  hardtopjet = PtEtaCut(300.0, 2.4); 
  topjetid = AndId<TopJet>(Type2TopTag(150,210, Type2TopTag::MassType::groomed,btag_medium),Tau32());
  ak4ForwardId = PtEtaCut(30.0,5,-1,2);
  ak4CentralId = PtEtaCut(30.0,2,-1,-1);

  common.reset(new CommonModules());
  common->set_jet_id(wide_softjet);
  common->set_electron_id(softElectron);
  common->set_muon_id(softMuon);
  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  common->init(ctx);
  common->set_HTjetid(softjet);

  channelSel.add<NElectronSelection>("0Electrons",0,0);
  channelSel.add<NMuonSelection>("1Muon",1,1,muid_cut);

  muonFactory.reset(new HistFactory(ctx));
  muonFactory->setEffiHistName("muonEffis");
  //muonFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*"),"muonTrigger");
  muonFactory->addSelection(make_unique<NElectronSelection>(0,0,softElectron),"0_eleCut");
  muonFactory->addSelection(make_unique<NMuonSelection>(1,1,softMuon),"1_softMuonCut");
  muonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  muonFactory->addSelection(make_unique<TwoDCut>(0.4,20),"2DCut");
  //muonFactory->addAnalysisModule(make_unique<JetCleaner>(secondjet));
  muonFactory->addSelection(make_unique<NMuonSelection>(1,1,muid_cut),"1_muonCut");
  muonFactory->addSelection(make_unique<NJetSelection>(2,-1,secondjet),"50GeV_JetCut");
  //muonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"130GeV_JetCut");
  muonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,hardtopjet),"300GeV_TopJetCut");
  //muonFactory->addSelection(make_unique<NTopJetSelection>(2,-1,topjet),"150GeV_2TopJetCut");
  //muonFactory->addSelection(make_unique<TopJetMassCut>(50),"50GeV_TopJetPruned");
  //muonFactory->addSelection(make_unique<ForwardJetPtEtaCut>(1.5,40),"ForwardJetCut");
  //muonFactory->addSelection(make_unique<STSelection>(ctx,800),"800GeV_ST");
  //muonFactory->addSelection(make_unique<NSubJetCut>(2),"2_SubJetCut");
  muonFactory->addSelection(make_unique<METSelection>(50),"50GeV_METCut");
  muonFactory->addSelection(make_unique<HTLepSelection>(ctx,150),"150_HTLep");
  //muonFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<HTLepSelection>(ctx,200)),"200_HTLep_or_TopTag");
  //muonFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<NJetSelection>(4,-1,softjet)),"4_JetCut_or_TopTag");
  muonFactory->addSelection(make_unique<NJetSelection>(4,-1,softjet),"4_JetCut");

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
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  //topWMuonFactory->addSelection(make_unique<TwoDCut>(0.4,40),"2DCut");
  //topWMuonFactory->addAnalysisModule(make_unique<JetCleaner>(secondjet));
  topWMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  topWMuonFactory->addSelection(make_unique<METSelection>(70),"70GeV_METCut");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,secondjet),"50GeV_JetCut");
  //topWMuonFactory->addSelection(make_unique<ForwardJetPtEtaCut>(1.5,40),"ForwardJetCut");
  //topWMuonFactory->addSelection(make_unique<STSelection>(ctx,900),"900GeV_ST");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"130GeV_JetCut");
  topWMuonFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<HTLepSelection>(ctx,200)),"200_HTLep_or_TopTag");
  topWMuonFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<NJetSelection>(4,-1,softjet)),"4_JetCut_or_TopTag");
  

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
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,softjet),"30GeV_JetCut");
  //wMuonFactory->addSelection(make_unique<TwoDCut>(0.4,40),"2DCut");
  //wMuonFactory->addAnalysisModule(make_unique<JetCleaner>(secondjet));
  wMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  wMuonFactory->addSelection(make_unique<METSelection>(70),"70GeV_METCut");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,secondjet),"50GeV_JetCut");
  //wMuonFactory->addSelection(make_unique<ForwardJetPtEtaCut>(1.5,40),"ForwardJetCut");
  wMuonFactory->addSelection(make_unique<STSelection>(ctx,900),"900GeV_ST");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"130GeV_JetCut");
  wMuonFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<HTLepSelection>(ctx,200)),"200_HTLep_or_TopTag");
  wMuonFactory->addOrSelection(make_uvec(make_unique<NTopJetSelection>(1,-1,topjetid),make_unique<NJetSelection>(4,-1,softjet)),"4_JetCut_or_TopTag");

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

bool GenTestModule::process(Event & event){
  if(!common->process(event)) return false;
  
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
 
  bool muonChannel = muonFactory->passAndFill(event); 
  return muonChannel;
 
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(GenTestModule)
