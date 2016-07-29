#include "UHH2/VLQToTopAndLepton/include/TopTagScalefactor.h"

using namespace uhh2;
using namespace std;


TopTagScalefactor::TopTagScalefactor(uhh2::Context & ctx, const std::string & hyp_name){
  recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
  weight_toptag = ctx.declare_event_output<double>("weight_toptag");
  weight_toptag_up = ctx.declare_event_output<double>("weight_toptag_up");
  weight_toptag_down = ctx.declare_event_output<double>("weight_toptag_down");
}


bool TopTagScalefactor::process(uhh2::Event & event){
  double weight = event.weight;
  double toptag_weight = -1;
  BprimeContainer hyp = event.get(recohyp);
  int recotype = hyp.get_RecoTyp();
  LorentzVector thad = hyp.get_topHad();
  
  if(recotype ==2 && !event.isRealData){
    if(thad.pt()>=400 && thad.pt()<=550){
      //event.weight *= 0.91;
      weight *= 0.91;
      toptag_weight = 0.91;
    }
    else if(thad.pt()>550){
      //event.weight *= 0.97;
      weight *= 0.97;
      toptag_weight = 0.97;
    }
  }

  event.set(weight_toptag,toptag_weight);
  if(toptag_weight>0){
    event.set(weight_toptag_up,toptag_weight*toptag_weight);
    event.set(weight_toptag_down,event.weight);
  }
  else{
    event.set(weight_toptag_up,-1);
    event.set(weight_toptag_down,-1);
  }
  return true;
}
