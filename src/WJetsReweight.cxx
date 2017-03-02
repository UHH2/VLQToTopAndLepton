#include "UHH2/VLQToTopAndLepton/include/WJetsReweight.h"


using namespace uhh2;
using namespace std;


WJetsReweight::WJetsReweight(Context & ctx, std::string mode_){
  string sample = ctx.get("dataset_version");
  mode = mode_;
  if(sample.find("WJet") != std::string::npos && sample.find("HT") != std::string::npos){
    work = true;
    wpt = new TF1("linear","pol1");
    h_wjets_reweight = ctx.declare_event_output<float>("weight_wjets_ht");
    h_wjets_reweight_up = ctx.declare_event_output<float>("weight_wjets_ht_up");
    h_wjets_reweight_down = ctx.declare_event_output<float>("weight_wjets_ht_down");
  }
  else
    work = false;
}


bool WJetsReweight::process(Event & event){
  if(!work || event.isRealData) return true;
  double wpt_val =0;
  for(auto genp : *event.genparticles){
    if(abs(genp.pdgId())==24){
      wpt_val = genp.pt();
      break;
    }
  }
  if(wpt_val<820){
    wpt->SetParameters(plin);
    if(mode=="up")
      wpt->SetParameters(plin_up);
    if(mode=="down")
      wpt->SetParameters(plin_down);
    event.weight *= wpt->Eval(wpt_val);

    wpt->SetParameters(plin);
    event.set(h_wjets_reweight, wpt->Eval(wpt_val));
    wpt->SetParameters(plin_up);
    event.set(h_wjets_reweight_up,0.889614);
    wpt->SetParameters(plin_down);
    event.set(h_wjets_reweight_down,0.564834);
  }
  else{
    if(mode=="central")
      event.weight *= 0.727;
    else if(mode=="up")
      event.weight *=0.889614;
    else if(mode=="down")
      event.weight *=0.564834;

    event.set(h_wjets_reweight,0.727);
    event.set(h_wjets_reweight_up,0.889614);
    event.set(h_wjets_reweight_down,0.564834);
  }
  return true;
}
