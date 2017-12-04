#include "UHH2/VLQToTopAndLepton/include/WtagFactors.h"



WtagFactors::WtagFactors(Context & ctx, TopJetId wtag_, string jetcollname){
  TFile * corrections = new TFile("weights/puppiCorr.root","READ");
  puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
  puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");
  corrections->Close();
  topjet_collection = ctx.get_handle<vector<TopJet>>(jetcollname);
  weight = ctx.declare_event_output<double>("weight_wtag");
}


bool WtagFactors::process(uhh2::Event & event){
  for(auto jet : event.get(topjet_collection)){
    
    float genCorr  = 1.;
    float recoCorr = 1.;
    float totalWeight = 1.;
    
    genCorr =  puppisd_corrGEN->Eval( jet.pt() );
    if( fabs(jet.eta())  <= 1.3 ){
      recoCorr = puppisd_corrRECO_cen->Eval( jet.pt() );
    }
    else{
      recoCorr = puppisd_corrRECO_for->Eval( jet.pt() );
    }
    for(auto subjet : jet.subjets()){
      subjet.v4().SetM(subjet.v4().M() * genCorr * recoCorr);
    }
  }
}
