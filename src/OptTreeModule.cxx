#include "UHH2/VLQToTopAndLepton/include/OptTreeModule.h"

using namespace uhh2; 
using namespace std;

OptTreeModule::OptTreeModule(Context & ctx){
  //Handles that need to be loaded 
  
  chi2HypHandle = ctx.get_handle<BprimeContainer>("Chi2Dis");
  TopTagHypHandle = ctx.get_handle<BprimeContainer>("TopTagDis");
  ttbarHandle = ctx.get_handle<BprimeContainer>("TTbarDis");
  //wTagHypHandle = ctx.get_handle<BprimeContainer>("WTagDis");
  //Output Variables 
  MET = ctx.declare_event_output<double>("MET");
  HTLep = ctx.declare_event_output<double>("HTLep");
  leadingLepPt = ctx.declare_event_output<double>("leadingLepPt");
  leadingJetPt = ctx.declare_event_output<double>("leadingJetPt");
  leadingTopJetPt = ctx.declare_event_output<double>("leadingTopJetPt");
  mostForwardJetEta = ctx.declare_event_output<double>("mostForwardJetEta");
  numberJet = ctx.declare_event_output<double>("numberJets");
  chi2Mass = ctx.declare_event_output<double>("chi2Mass");
  chi2Val = ctx.declare_event_output<double>("chi2Val");
  TopTagMass = ctx.declare_event_output<double>("TopTagMass");
  TopTagChi2 = ctx.declare_event_output<double>("TopTagChi2");
  ttbarMass = ctx.declare_event_output<double>("ttbarMass");
  ttbarChi2 = ctx.declare_event_output<double>("ttbarChi2");
  ttbarMass = ctx.declare_event_output<double>("ttbarMass");
  weight = ctx.declare_event_output<double>("weight");
  pTW = ctx.declare_event_output<double>("pTW");
  pTT = ctx.declare_event_output<double>("pTT");
  NBmediumtags = ctx.declare_event_output<double>("NBmediumtags");
  NBtighttags = ctx.declare_event_output<double>("NBtighttags");
  NBloosetags = ctx.declare_event_output<double>("NBloosetags");
  mediumbtag_id = CSVBTag(CSVBTag::WP_MEDIUM);
  tightbtag_id = CSVBTag(CSVBTag::WP_TIGHT);
  loosebtag_id = CSVBTag(CSVBTag::WP_LOOSE);
  
}

bool OptTreeModule::process(Event & event){
  event.set(MET,event.met->pt());
  event.set(leadingJetPt,event.jets->at(0).pt());
  event.set(leadingTopJetPt,event.jets->at(0).pt());
  if(event.muons->size()>0){
    event.set(HTLep,event.met->pt()+event.muons->at(0).pt());
    event.set(leadingLepPt,event.muons->at(0).pt());
  }
  else{
    event.set(HTLep,event.met->pt()+event.electrons->at(0).pt());
    event.set(leadingLepPt,event.electrons->at(0).pt());
  }
  event.set(numberJet,event.jets->size());
  LorentzVector forwardjet(0,0,0,0);
  int n_btagmedium =0;
  int n_btagtight=0;
  int n_btagloose=0;
  for(auto jet : *event.jets){
    if((mediumbtag_id)(jet, event)) ++n_btagmedium;
    if((tightbtag_id)(jet, event)) ++n_btagtight;
    if((loosebtag_id)(jet, event)) ++n_btagloose;
    if(jet.eta()>forwardjet.eta())
      forwardjet = jet.v4();
  }
  event.set(NBmediumtags,n_btagmedium);
  event.set(NBtighttags,n_btagtight);
  event.set(NBloosetags,n_btagloose);
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
  event.set(TopTagMass,event.get(TopTagHypHandle).get_Mass());
  event.set(TopTagChi2,event.get(TopTagHypHandle).get_chiVal());
  //event.set(wTagMass,event.get(wTagHypHandle).get_Mass());
  return true;
}
