#include "UHH2/VLQToTopAndLepton/include/OptTreeModule.h"

using namespace uhh2;
using namespace std;

OptTreeModule::OptTreeModule(Context & ctx){
  leadingJetPt = ctx.declare_event_output<double>("leadingJetPt");
  numberJet = ctx.declare_event_output<double>("numberJet");
  /*
  chi2Mass = ctx.declare_event_output<double>("chi2Mass");
  chi2Val = ctx.declare_event_output<double>("chi2Val");
  ttbarchi2 = ctx.declare_event_output<double>("ttbarChi2");
  cmsTopTagMass = ctx.declare_event_output<double>("cmsTopTagMass");
  cmsTopTagChi2 = ctx.declare_event_output<double>("cmsTopTagChi2");
  pTW = ctx.declare_event_output<double>("pTW");
  pTT = ctx.declare_event_output<double>("pTT");
  */

}

bool OptTreeModule::process(Event & event){
  event.set(leadingJetPt,event.jets->at(0).pt());
  event.set(numberJet,event.jets->size());
  //event.set();
  //event.set();
  //event.set();
  //event.set();
  //event.set();
  return true;
}
