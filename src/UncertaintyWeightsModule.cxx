#include "UHH2/VLQToTopAndLepton/include/UncertaintyWeightsModule.h"

using namespace std;
using namespace uhh2;

UncertaintyWeightsModule::UncertaintyWeightsModule(uhh2::Context & ctx){
  scaleWeight_up = ctx.declare_event_output<double>("scaleWeight_up");
  scaleWeight_down = ctx.declare_event_output<double>("scaleWeight_down");
  pdfWeight = ctx.declare_event_output<double>("pdfWeight");

}

bool UncertaintyWeightsModule::process(uhh2::Event & event){
  double pdf_final = 0;
  if(event.isRealData){
    event.set(scaleWeight_up,-1);
    event.set(scaleWeight_down,-1);
    event.set(pdfWeight,-1);
    return true;
  }

  double pdf_weight = event.genInfo->weights().at(0);
  for(unsigned int i = 9; i<109; i++){
    if(event.genInfo->systweights().size()==0){ 
      pdf_final = -1;
      break;
    }
    double tmp_weight = event.genInfo->systweights().at(i)/event.genInfo->originalXWGTUP();
    pdf_final += (pdf_weight-tmp_weight)*(pdf_weight-tmp_weight);
  }
  pdf_final = sqrt(pdf_final/99);

  double scale_final_up = -9999999999;
  double scale_final_down = 999999999;
  for(auto entry : scale_entries){
    if(event.genInfo->systweights().size()==0){
      scale_final_up =-1;
      scale_final_down =-1;
      break;
    }
    double tmp_weight = event.genInfo->systweights().at(entry)/event.genInfo->originalXWGTUP();
    if(tmp_weight > scale_final_up)
      scale_final_up = tmp_weight;
    if(tmp_weight < scale_final_down)
      scale_final_down = tmp_weight;
  }
  //cout<<"nominal "<<event.weight<<" scale up "<<scale_final_up*event.weight<<" scale down "<<scale_final_down*event.weight<<endl;

  event.set(scaleWeight_up,scale_final_up);
  event.set(scaleWeight_down,scale_final_down);
  event.set(pdfWeight,pdf_final);
  return true;
}
