#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Selection.h"

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

using namespace std;
using namespace uhh2;


class GenTestModule: public AnalysisModule {
public:

  explicit GenTestModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;

  std::unique_ptr<VLQGenHists> vlqGenHists;
  std::unique_ptr<AnalysisModule> ht;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> bBprimeFactory;
  std::unique_ptr<HistFactory> muonFactory, dileptonFactory;//, hadronicFactory;
  std::unique_ptr<HistFactory> topWMuonFactory, wMuonFactory;
  std::unique_ptr<HistFactory> muonTrigger;
  std::unique_ptr<AnalysisModule>  jetCorr; 
  std::unique_ptr<ElectronCleaner> elecleaner;
  std::unique_ptr<MuonCleaner> mucleaner;  
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<AnalysisModule> jetlepclean; 
  std::unique_ptr<Selection> HiggsFilter, ZFilter;
  AndSelection  channelSel;
  JetId btag_medium;
  MuonId muid_cut;
  JetId twojet, onejet, lowjet;
  TopJetId topjet;
  //std::unique_ptr<BprimeReco> Reco;
};



GenTestModule::GenTestModule(Context& ctx):channelSel(ctx){

  //Reco.reset(new BprimeReco(ctx)); 
  Version  = ctx.get("dataset_version", "<not set>");
  btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  ht.reset(new HTCalc(ctx));
  lepton.reset(new PrimaryLepton(ctx));

  jetcleaner.reset(new JetCleaner(30.0, 3));
  jetlepclean.reset(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));
  jetCorr.reset(new JetCorrector(JERFiles::PHYS14_L123_MC));
  elecleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_CSA14_50ns_medium, PtEtaCut(35.0, 2.4))));
  mucleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4))));

  vlqGenHists.reset(new VLQGenHists(ctx,"GenHists"));

  muonTrigger.reset(new HistFactory(ctx,"triggerEffis.txt"));
  //muonTrigger.reset(new HistFactory(ctx));
  muonTrigger->addSelection(make_unique<GenNSelection>(13,1,1,30,-1),"1_GenSel");
  muonTrigger->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"muonTrigger");
  muonTrigger->addHists("MuonHists","triggerChannel_MuonHists");

  muid_cut = AndId<Muon>(MuonIDTight(), PtEtaCut(50.0, 2.4));
  onejet = PtEtaCut(200.0, 2.4); twojet = PtEtaCut(50.0, 2.4);
  topjet = PtEtaCut(250.0,2.4); 
  lowjet = PtEtaCut(50.0,2.4); 

  channelSel.add<NElectronSelection>("0Electrons",0,0);
  channelSel.add<NMuonSelection>("1Muon",1,1,muid_cut);


  //muonFactory.reset(new HistFactory(ctx,"muonEffis.txt"));
  muonFactory.reset(new HistFactory(ctx));
  muonFactory->setEffiHistName("muonEffis");
  //muonFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"muonTrigger");
  muonFactory->addSelection(make_unique<NElectronSelection>(0,0),"0_eleCut");
  muonFactory->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  //muonFactory->addSelection(make_unique<STSelection>(ctx,500),"500_ST");
  muonFactory->addSelection(make_unique<NJetSelection>(1,-1,lowjet),"50GeV_JetCut");
  muonFactory->addSelection(make_unique<TwoDCut>(0.4,25),"2DCut");
  muonFactory->addSelection(make_unique<HTLepSelection>(ctx,200),"200_HTLep");
  //muonFactory->addSelection(make_unique<METSelection>(ctx,50),"50_MET");
  muonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"200GeV_JetCut");
  muonFactory->addSelection(make_unique<NJetSelection>(2,-1,twojet),"3_JetCut");
  //muonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,topjet), "TopJetCut");
  //muonFactory->addSelection(make_unique<NJetSelection>(1,-1, btag_medium), "BTagMedium");

  muonFactory->addHists("ElectronHists","muonChannel_ElectronHists");
  muonFactory->addHists("MuonHists","muonChannel_MuonHists");
  muonFactory->addHists("EventHists","muonChannel_EventHists");
  muonFactory->addHists("JetHists","muonChannel_JetHists");
  muonFactory->addHists("TopJetHists","muonChannel_TopJetHists");
  muonFactory->addHists("VLQGenHists","muonChannel_VLQGenHists");


  /*
  dileptonFactory.reset(new HistFactory(ctx,"dileptonEffis.txt"));	
  //dileptonFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"muonTrigger");
  dileptonFactory->addSelection(make_unique<NElectronSelection>(0,0),"0_eleCut");
  dileptonFactory->addSelection(make_unique<NMuonSelection>(2,2,muid_cut),"2_muonCut");
  dileptonFactory->addSelection(make_unique<NJetSelection>(2),"2_JetCut");
  dileptonFactory->addSelection(make_unique<NJetSelection>(2,-1, btag_medium), "BTagMedium");

  dileptonFactory->addHists("ElectronHists","dileptonChannel_ElectronHists");
  dileptonFactory->addHists("MuonHists","dileptonChannel_MuonHists");
  dileptonFactory->addHists("EventHists","dileptonChannel_EventHists");
  dileptonFactory->addHists("JetHists","dileptonChannel_JetHists");
  dileptonFactory->addHists("TopJetHists","dileptonChannel_TopJetHists");
  dileptonFactory->addHists("VLQGenHists","dileptonChannel_VLQGenHists");
  */
  /*
  hadronicFactory.reset(new HistFactory(ctx,"hadronicEffis.txt"));	
  hadronicFactory->addSelection(make_unique<NMuonSelection>(0,0),"0_muonCut");
  hadronicFactory->addSelection(make_unique<NElectronSelection>(0,0),"0_eleCut");
  hadronicFactory->addSelection(make_unique<NJetSelection>(2),"2_JetCut");
  hadronicFactory->addSelection(make_unique<NTopJetSelection>(2),"2_TopJetCut");
  hadronicFactory->addSelection(make_unique<NJetSelection>(2,-1, btag_medium), "BTagMedium");

  hadronicFactory->addHists("ElectronHists","hadronicChannel_ElectronHists");
  hadronicFactory->addHists("MuonHists","hadronicChannel_MuonHists");
  hadronicFactory->addHists("EventHists","hadronicChannel_EventHists");
  hadronicFactory->addHists("JetHists","hadronicChannel_JetHists");
  hadronicFactory->addHists("TopJetHists","hadronicChannel_TopJetHists");
  hadronicFactory->addHists("VLQGenHists","hadronicChannel_VLQGenHists");
  */

  vector<int> bBprimeDecay {5, 6000007};

  bBprimeFactory.reset(new HistFactory(ctx,"bBprime.txt"));
  bBprimeFactory->addSelection(make_unique<GenFamilySelection>(bBprimeDecay,2),"bB'");


  vector<int> topLep {6,24,13};
  vector<int> topHad {6,24,-54321};
  vector<int> wHad {24,-54321};
  vector<int> wLep {24,13};


  //topWMuonFactory.reset(new HistFactory(ctx,"topLep.txt"));
  topWMuonFactory.reset(new HistFactory(ctx));
  topWMuonFactory->setEffiHistName("topLep");
  topWMuonFactory->addSelection(make_unique<GenFamilySelection>(topLep,2),"topLep");
  topWMuonFactory->addSelection(make_unique<GenFamilySelection>(wHad,2),"wHad");
  topWMuonFactory->addSelection(make_unique<NElectronSelection>(0,0),"0_eleCut");
  topWMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  //topWMuonFactory->addSelection(make_unique<STSelection>(ctx,500),"500_ST");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,lowjet),"50GeV_JetCut");
  topWMuonFactory->addSelection(make_unique<TwoDCut>(0.4,25),"2DCut");
  topWMuonFactory->addSelection(make_unique<HTLepSelection>(ctx,150),"200_HTLep");
  //topWMuonFactory->addSelection(make_unique<METSelection>(ctx,50),"50_MET");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"200GeV_JetCut");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(2,-1,twojet),"3_JetCut");
  //topWMuonFactory->addSelection(make_unique<NTopJetSelection>(1,-1,topjet), "TopJetCut");
  topWMuonFactory->addSelection(make_unique<NJetSelection>(1,-1, btag_medium), "BTagMedium");

  topWMuonFactory->addHists("ElectronHists","topLep_ElectronHists");
  topWMuonFactory->addHists("MuonHists","topLep_MuonHists");
  topWMuonFactory->addHists("EventHists","topLep_EventHists");
  topWMuonFactory->addHists("JetHists","topLep_JetHists");
  topWMuonFactory->addHists("TopJetHists","topLep_TopJetHists");
  topWMuonFactory->addHists("VLQGenHists","topLep_VLQGenHists");

  //wMuonFactory.reset(new HistFactory(ctx,"topHad.txt"));
  wMuonFactory.reset(new HistFactory(ctx));
  wMuonFactory->setEffiHistName("topHad");	
  wMuonFactory->addSelection(make_unique<GenFamilySelection>(topHad,2),"topHad");
  wMuonFactory->addSelection(make_unique<GenFamilySelection>(wLep,2),"wLep");
  wMuonFactory->addSelection(make_unique<NElectronSelection>(0,0),"0_eleCut");
  wMuonFactory->addSelection(make_unique<NMuonSelection>(1,-1,muid_cut),"1_muonCut");
  //wMuonFactory->addSelection(make_unique<STSelection>(ctx,500),"500_ST");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,lowjet),"50GeV_JetCut");
  wMuonFactory->addSelection(make_unique<TwoDCut>(0.4,25),"2DCut");
  wMuonFactory->addSelection(make_unique<HTLepSelection>(ctx,200),"200_HTLep");
  //wMuonFactory->addSelection(make_unique<METSelection>(ctx,50),"50_MET");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1,onejet),"200GeV_JetCut");
  wMuonFactory->addSelection(make_unique<NJetSelection>(2,-1,twojet),"3_JetCut");
  //wMuonFactory->addSelection(make_unique<NTopJetSelection>(2,-1,topjet), "TopJetCut");
  wMuonFactory->addSelection(make_unique<NJetSelection>(1,-1, btag_medium), "BTagMedium");

  wMuonFactory->addHists("ElectronHists","topHad_ElectronHists");
  wMuonFactory->addHists("MuonHists","topHad_MuonHists");
  wMuonFactory->addHists("EventHists","topHad_EventHists");
  wMuonFactory->addHists("JetHists","topHad_JetHists");
  wMuonFactory->addHists("TopJetHists","topHad_TopJetHists");
  wMuonFactory->addHists("VLQGenHists","topHad_VLQGenHists");


  HiggsFilter.reset(new GenParticleFilter(25,0,0));
  ZFilter.reset(new GenParticleFilter(23,0,0));

}

bool GenTestModule::process(Event & event){
  elecleaner->process(event);
  mucleaner->process(event);
  jetlepclean->process(event);
  jetcleaner->process(event);
  jetCorr->process(event);
  ht->process(event);
  lepton->process(event);
  sort_by_pt(*event.jets);

  bBprimeFactory->passAndFill(event);
  vlqGenHists->fill(event);    
 
  if(Version.find("BpJ_TW") != std::string::npos)
    if(!HiggsFilter->passes(event)||!ZFilter->passes(event)) return false;

  //Reco->massReco(event);
  //channelSel.passes(event);  
  //muonTrigger->passAndFill(event);
  wMuonFactory->passAndFill(event);
  topWMuonFactory->passAndFill(event);

  //hadronicFactory->passAndFill(event);
  //dileptonFactory->passAndFill(event);
  bool muonChannel = muonFactory->passAndFill(event); 
  return muonChannel;
 
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(GenTestModule)
