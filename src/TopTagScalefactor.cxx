#include "UHH2/VLQToTopAndLepton/include/TopTagScalefactor.h"

using namespace uhh2;
using namespace std;


TopTagScalefactor::TopTagScalefactor(uhh2::Context & ctx, const std::string & hyp_name, const std::string & wtag ){
  wtaghyp = ctx.get_handle<BprimeContainer>(wtag);   
  recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
  
  weight_toptag = ctx.declare_event_output<double>("weight_toptag");
  weight_toptag_up = ctx.declare_event_output<double>("weight_toptag_up");
  weight_toptag_down = ctx.declare_event_output<double>("weight_toptag_down");
  
  weight_wtag = ctx.declare_event_output<double>("weight_wtag");
  weight_wtag_up = ctx.declare_event_output<double>("weight_wtag_up");
  weight_wtag_down = ctx.declare_event_output<double>("weight_wtag_down");
}


bool TopTagScalefactor::process(uhh2::Event & event){
  if(event.isRealData){
    event.set(weight_wtag,event.weight);
    event.set(weight_wtag_up,event.weight);
    event.set(weight_wtag_down,event.weight);
  }
  BprimeContainer whyp = event.get(wtaghyp);
  if(whyp.get_wHad().pt()>200){
    event.set(weight_wtag,event.weight*0.91);
    event.set(weight_wtag_up,event.weight*0.96);
    event.set(weight_wtag_down,event.weight*.86);
  }
  else{
    event.set(weight_wtag,-1.);
    event.set(weight_wtag_up,-1.);
    event.set(weight_wtag_down,-1.);
  }
  BprimeContainer tophyp = event.get(recohyp);
  if(tophyp.get_topHad().pt()>100){
    event.set(weight_toptag,event.weight*1.01);
    event.set(weight_toptag_up,event.weight*1.08);
    event.set(weight_toptag_down,event.weight*.97);
  }
  else{
    event.set(weight_toptag,-1.);
    event.set(weight_toptag_up,-1.);
    event.set(weight_toptag_down,-1.);
  }
 
  return true;
}
