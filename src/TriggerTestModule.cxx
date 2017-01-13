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

using namespace std;
using namespace uhh2;


class TriggerTestModule: public AnalysisModule {
public:

  explicit TriggerTestModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  string Version;

  std::unique_ptr<AnalysisModule> ht;
  std::unique_ptr<AnalysisModule> lepton;
  std::unique_ptr<HistFactory> muonTriggerFactory, muonGenTriggerFactory;
  std::unique_ptr<HistFactory> muonTriggerIsoFactory, muonGenTriggerIsoFactory;
  std::unique_ptr<HistFactory> muonTriggerIsoSelFactory, muonTriggerSelFactory;
  std::unique_ptr<AnalysisModule>  jetCorr; 
  std::unique_ptr<ElectronCleaner> elecleaner;
  std::unique_ptr<MuonCleaner> mucleaner, muIsoCleaner;  
  std::unique_ptr<JetCleaner> jetcleaner;
  std::unique_ptr<Selection> recoMuon, genMuon;
  std::unique_ptr<Selection> HiggsFilter, ZFilter;
  std::unique_ptr<Selection> ElectronSelection, MuonSelection, IsoMuonSelection, stSelection, JetSelection_one, JetSelection_two;
  MuonId muid_cut, isomuon;
  JetId twojet, onejet;
};



TriggerTestModule::TriggerTestModule(Context& ctx){


  Version  = ctx.get("dataset_version", "<not set>");

  ht.reset(new HTCalc(ctx));
  lepton.reset(new PrimaryLepton(ctx));

  genMuon.reset(new GenNSelection(13,1,1,20,-1));
  recoMuon.reset(new NMuonSelection(1,1));

  jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4));
  //jetCorr.reset(new JetCorrector(ctx,JERFiles::Summer15_25ns_L123_AK4PFchs_MC));
  elecleaner.reset(new ElectronCleaner(AndId<Electron>(ElectronID_Spring15_50ns_medium, PtEtaCut(30.0, 2.4))));
  mucleaner.reset(new MuonCleaner(AndId<Muon>(MuonIDTight(),PtEtaCut(20.0, 3.0))));

  muIsoCleaner.reset(new MuonCleaner(MuonIso()));


  muid_cut = AndId<Muon>(MuonIDTight(), PtEtaCut(50.0, 2.4)); onejet = PtEtaCut(250.0, 2.4); twojet = PtEtaCut(50.0, 2.4);isomuon = MuonIso();
  ElectronSelection.reset(new NElectronSelection(0,0)); MuonSelection.reset(new NMuonSelection(1,1,muid_cut));stSelection.reset(new STSelection(ctx,500));JetSelection_one.reset(new NJetSelection(1,-1,onejet));JetSelection_two.reset(new NJetSelection(3,-1,twojet));IsoMuonSelection.reset(new NMuonSelection(1,1,isomuon));



  muonGenTriggerFactory.reset(new HistFactory(ctx,"GenTriggerEffis.txt"));
  //muonGenTriggerFactory.reset(new HistFactory(ctx));
  muonGenTriggerFactory->setEffiHistName("muonGenTriggerEffis");
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v1"),"Mu40");
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v1"),"Mu40Jets");
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35");
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40");		   
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1"),"IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV");  
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_v1"),"IsoMu24_eta2p1_IterTrk02");						   
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_eta2p1_IterTrk02_v1"),"IsoTkMu24_eta2p1_IterTrk02");						   
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_IterTrk02_v1"),"IsoMu24_IterTrk02");								   
  muonGenTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_IterTrk02_v1"),"IsoTkMu24_IterTrk02");                                                                

  muonGenTriggerFactory->addHists("MuonHists","MuonGen_MuonHists");
  muonGenTriggerFactory->addHists("JetHists","MuonGen_JetHists");
  muonGenTriggerFactory->addHists("VLQGenHists","MuonGen_VLQGenHists");     

   
  muonTriggerFactory.reset(new HistFactory(ctx,"RecoTriggerEffis.txt"));
  //muonTriggerFactory.reset(new HistFactory(ctx));
  muonTriggerFactory->setEffiHistName("muonTriggerEffis");
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"Mu40");
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v1"),"Mu40Jets");
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35");
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40");		   
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1"),"IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV");  
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_v1"),"IsoMu24_eta2p1_IterTrk02");						     
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_eta2p1_IterTrk02_v1"),"IsoTkMu24_eta2p1_IterTrk02");						     
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_IterTrk02_v1"),"IsoMu24_IterTrk02");								   
  muonTriggerFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_IterTrk02_v1"),"IsoTkMu24_IterTrk02");          


  muonTriggerFactory->addHists("MuonHists","MuonReco_MuonHists");
  muonTriggerFactory->addHists("JetHists","MuonReco_JetHists");
  muonTriggerFactory->addHists("VLQGenHists","MuonReco_VLQGenHists");   
    
 

  muonGenTriggerIsoFactory.reset(new HistFactory(ctx,"GenTriggerEffisIso.txt"));
  //muonGenTriggerIsoFactory.reset(new HistFactory(ctx));
  muonGenTriggerIsoFactory->setEffiHistName("muonGenTriggerEffisIso");
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v1"),"Mu40");
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v1"),"Mu40Jets");
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35");
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40");		   
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1"),"IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV");  
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_v1"),"IsoMu24_eta2p1_IterTrk02");						   
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_eta2p1_IterTrk02_v1"),"IsoTkMu24_eta2p1_IterTrk02");						   
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_IterTrk02_v1"),"IsoMu24_IterTrk02");								   
  muonGenTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_IterTrk02_v1"),"IsoTkMu24_IterTrk02");                                                                

  muonGenTriggerIsoFactory->addHists("MuonHists","MuonGenIso_MuonHists");
  muonGenTriggerIsoFactory->addHists("JetHists","MuonGenIso_JetHists");
  muonGenTriggerIsoFactory->addHists("VLQGenHists","MuonGenIso_VLQGenHists");     

   
  muonTriggerIsoFactory.reset(new HistFactory(ctx,"RecoTriggerEffisIso.txt"));
  //muonTriggerIsoFactory.reset(new HistFactory(ctx));
  muonTriggerIsoFactory->setEffiHistName("muonTriggerEffisIso");
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"Mu40");
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v1"),"Mu40Jets");
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35");
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40");		   
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1"),"IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV");  
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_v1"),"IsoMu24_eta2p1_IterTrk02");						     
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_eta2p1_IterTrk02_v1"),"IsoTkMu24_eta2p1_IterTrk02");						     
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_IterTrk02_v1"),"IsoMu24_IterTrk02");								   
  muonTriggerIsoFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_IterTrk02_v1"),"IsoTkMu24_IterTrk02");          

  muonTriggerIsoFactory->addHists("MuonHists","MuonRecoIso_MuonHists");
  muonTriggerIsoFactory->addHists("JetHists","MuonRecoIso_JetHists");
  muonTriggerIsoFactory->addHists("VLQGenHists","MuonRecoIso_VLQGenHists");      


  muonTriggerSelFactory.reset(new HistFactory(ctx,"RecoTriggerSelEffis.txt"));
  //muonTriggerSelFactory.reset(new HistFactory(ctx));
  muonTriggerSelFactory->setEffiHistName("muonTriggerSelEffis");
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"Mu40");
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v1"),"Mu40Jets");
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35");
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40");		   
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1"),"IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV");  
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_v1"),"IsoMu24_eta2p1_IterTrk02");						     
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_eta2p1_IterTrk02_v1"),"IsoTkMu24_eta2p1_IterTrk02");						     
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_IterTrk02_v1"),"IsoMu24_IterTrk02");								   
  muonTriggerSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_IterTrk02_v1"),"IsoTkMu24_IterTrk02");          


  muonTriggerSelFactory->addHists("MuonHists","MuonRecoSel_MuonHists");
  muonTriggerSelFactory->addHists("JetHists","MuonRecoSel_JetHists");
  muonTriggerSelFactory->addHists("VLQGenHists","MuonRecoSel_VLQGenHists");   
  

  muonTriggerIsoSelFactory.reset(new HistFactory(ctx,"RecoTriggerEffisIsoSel.txt"));
  //muonTriggerIsoSelFactory.reset(new HistFactory(ctx));
  muonTriggerIsoSelFactory->setEffiHistName("muonTriggerEffisIsoSel");
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_v*"),"Mu40");
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_Mu40_eta2p1_PFJet200_PFJet50_v1"),"Mu40Jets");
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35");
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v1"),"IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40");		   
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1"),"IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV");  
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_eta2p1_IterTrk02_v1"),"IsoMu24_eta2p1_IterTrk02");						     
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_eta2p1_IterTrk02_v1"),"IsoTkMu24_eta2p1_IterTrk02");						     
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoMu24_IterTrk02_v1"),"IsoMu24_IterTrk02");								   
  muonTriggerIsoSelFactory->addSelection(make_unique<TriggerSelection>("HLT_IsoTkMu24_IterTrk02_v1"),"IsoTkMu24_IterTrk02");          


  muonTriggerIsoSelFactory->addHists("MuonHists","MuonRecoIsoSel_MuonHists");
  muonTriggerIsoSelFactory->addHists("JetHists","MuonRecoIsoSel_JetHists");
  muonTriggerIsoSelFactory->addHists("VLQGenHists","MuonRecoIsoSel_VLQGenHists");    


                                
  /*
    muonTriggerFactory->addHists("ElectronHists","muonChannel_ElectronHists");
    muonTriggerFactory->addHists("EventHists","muonChannel_EventHists");
    muonTriggerFactory->addHists("TopJetHists","muonChannel_TopJetHists");
   
  */

  HiggsFilter.reset(new GenParticleFilter(25,0,0));
  ZFilter.reset(new GenParticleFilter(23,0,0));

}

bool TriggerTestModule::process(Event & event){

  //jetCorr->process(event);
  jetcleaner->process(event);
  elecleaner->process(event);
  mucleaner->process(event);
 

  ht->process(event);
  lepton->process(event);
  
  
  if(Version.find("BpJ_TW") != std::string::npos)
    if(!HiggsFilter->passes(event)||!ZFilter->passes(event)) return false;
  
  
  if(ElectronSelection->passes(event) && MuonSelection->passes(event) && stSelection->passes(event) && JetSelection_one->passes(event) && JetSelection_two->passes(event)){
    muonTriggerSelFactory->passAndFill(event,1);
    if(IsoMuonSelection->passes(event)) muonTriggerIsoSelFactory->passAndFill(event,1);
  }

  if(genMuon->passes(event)) muonGenTriggerFactory->passAndFill(event,1);
  if(recoMuon->passes(event)) {
    muonTriggerFactory->passAndFill(event,1);   
  }
  else return false;

  muIsoCleaner->process(event);
  if(genMuon->passes(event))muonGenTriggerIsoFactory->passAndFill(event,1);
  if(recoMuon->passes(event)) muonTriggerIsoFactory->passAndFill(event,1);
  




  return true;
}


// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TriggerTestModule)
