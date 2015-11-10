#include "UHH2/VLQToTopAndLepton/include/OptTreeModule.h"

using namespace uhh2;
using namespace std;

OptTreeModule::OptTreeModule(Context & ctx){
  //Handles that need to be loaded 
  chi2HypHandle = ctx.get_handle<BprimeContainer>("Chi2Dis");
  cmsTopTagHypHandle = ctx.get_handle<BprimeContainer>("CMSTopTagDis");
  ttbarHandle = ctx.get_handle<BprimeContainer>("TTbarDis");
  wTagHypHandle = ctx.get_handle<BprimeContainer>("WTagDis");
  //Output Variables  
  leadingLepPt = ctx.declare_event_output<double>("leadingLepPt");
  leadingJetPt = ctx.declare_event_output<double>("leadingJetPt");
  mostForwardJetEta = ctx.declare_event_output<double>("mostForwardJetEta");
  numberJet = ctx.declare_event_output<double>("numberJets");
  chi2Mass = ctx.declare_event_output<double>("chi2Mass");
  chi2Val = ctx.declare_event_output<double>("chi2Val");
  cmsTopTagMass = ctx.declare_event_output<double>("cmsTopTagMass");
  cmsTopTagChi2 = ctx.declare_event_output<double>("cmsTopTagChi2");
  ttbarMass = ctx.declare_event_output<double>("ttbarMass");
  ttbarChi2 = ctx.declare_event_output<double>("ttbarChi2");
  ttbarMass = ctx.declare_event_output<double>("ttbarMass");
  weight = ctx.declare_event_output<double>("weight");
  pTW = ctx.declare_event_output<double>("pTW");
  pTT = ctx.declare_event_output<double>("pTT");
}

bool OptTreeModule::process(Event & event){
  event.set(leadingJetPt,event.jets->at(0).pt());
  event.set(numberJet,event.jets->size());
  LorentzVector forwardjet(0,0,0,0);
  for(auto jet : *event.jets){
    if(jet.eta()>forwardjet.eta())
      forwardjet = jet.v4();
  }
  event.set(mostForwardJetEta,forwardjet.eta());
  event.set(weight,event.weight);
  event.set(chi2Mass,event.get(chi2HypHandle).get_Mass());
  event.set(chi2Val,event.get(chi2HypHandle).get_chiVal());
  if(event.get(chi2HypHandle).get_Mass()!=-1){  
    if((event.get(chi2HypHandle)).get_RecoTyp() == 11)
      event.set(pTT,(event.get(chi2HypHandle)).get_topLep().pt()/(event.get(chi2HypHandle)).get_wHad().pt());
    else
      event.set(pTT,(event.get(chi2HypHandle)).get_topHad().pt()/(event.get(chi2HypHandle)).get_wLep().pt());
    event.set(pTW,(event.get(chi2HypHandle)).get_wLep().pt()/(event.get(chi2HypHandle)).get_wHad().pt());
  }
  else{
    event.set(pTW,-1);
    event.set(pTT,-1);
  }
  event.set(ttbarMass,event.get(ttbarHandle).get_Mass());
  event.set(ttbarChi2,event.get(ttbarHandle).get_chiVal());
  event.set(cmsTopTagMass,event.get(cmsTopTagHypHandle).get_Mass());
  event.set(cmsTopTagChi2,event.get(cmsTopTagHypHandle).get_chiVal());
  event.set(wTagMass,event.get(wTagHypHandle).get_Mass());
  return true;
}
