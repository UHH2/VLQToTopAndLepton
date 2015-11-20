#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"

#include "UHH2/common/include/Utils.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;


GenParticleFilter::GenParticleFilter(int pdgId_, int nmin_, int nmax_):pdgId(pdgId_), nmin(nmin_),nmax(nmax_)  {}

bool GenParticleFilter::passes(const uhh2::Event & event){  
  int count = 0 ;
  for(auto genp : *event.genparticles){
    if(pdgId == abs(genp.pdgId())) count++;
    if(count>nmax) return false;
  }
  return count<=nmax && count>=nmin;
}


HTSelection::HTSelection(Context & ctx, double HTmin_): HTmin(HTmin_){
  ht = ctx.get_handle<double>("HT");
}

bool HTSelection::passes(const Event & event){
  return HTmin<event.get(ht);
  //return false;
}

STSelection::STSelection(Context & ctx, double STmin_): STmin(STmin_){
  ht = ctx.get_handle<double>("HT");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
}

bool STSelection::passes(const Event & event){
  return STmin<(event.get(ht)+event.met->pt()+event.get(h_primlep).pt());
}

HTLepSelection::HTLepSelection(Context & ctx, double HTLepmin_): HTLepmin(HTLepmin_){
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
}

bool HTLepSelection::passes(const Event & event){
  return HTLepmin < event.met->pt()+event.get(h_primlep).pt();
}

METSelection::METSelection(double METmin_): METmin(METmin_){}

bool METSelection::passes(const Event & event){
  return METmin < event.met->pt();
}

bool TwoDCut::passes(const Event & event){

  assert((event.muons || event.electrons) && event.jets);/*
  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- TwoDCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    cout<<"muons "<<event.muons->size()<<" ele "<<event.electrons->size()<<endl;
    for(auto muon : *event.muons)
      cout<<muon.pt()<<" ";
    cout<<endl;
    return false;
  }
						       */
  float drmin, ptrel;  
  if(event.muons->size()) std::tie(drmin, ptrel) = drmin_pTrel(event.muons->at(0), *event.jets);
  else std::tie(drmin, ptrel) = drmin_pTrel(event.electrons->at(0), *event.jets);
  //cout<<"dR "<<drmin<<" pTrel "<<ptrel<<endl;

  return (drmin > min_deltaR_) || (ptrel > min_pTrel_);
}


ChiSquareCut::ChiSquareCut(Context & ctx, float max_chi2,float min_chi2, const std::string & hyp_name, int recotyp):min_(min_chi2), max_(max_chi2),recotyp_(recotyp){
  recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
}

bool ChiSquareCut::passes(const uhh2::Event & event){
  BprimeContainer hyp = event.get(recohyp);
  if(recotyp_==-1 || recotyp_ ==  hyp.get_RecoTyp()){
    if((max_>hyp.get_chiVal() && hyp.get_chiVal()>min_) || (max_>hyp.get_chiVal() && min_==-1) || (min_<hyp.get_chiVal() && max_==-1)) return true;
  }
  else if(recotyp_!= hyp.get_RecoTyp()) return true;
  return false;
}

PtRatioWTCut::PtRatioWTCut(Context & ctx, float min, float max, const std::string & hyp_name):min_(min), max_(max){
  recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
}

bool PtRatioWTCut::passes(const uhh2::Event & event){
  BprimeContainer hyp = event.get(recohyp);
  float ratio =-1;
  if(hyp.get_RecoTyp()==11){
    ratio = hyp.get_wHad().pt()/hyp.get_topLep().pt();
  }
  else if(hyp.get_RecoTyp()==12 || hyp.get_RecoTyp()==2){
    ratio = hyp.get_wLep().pt()/hyp.get_topHad().pt();
  }
  if(ratio==-1){
    cout<<" pT ratio W/Top is -1. Have a look at the cut?"<<endl;
    exit(EXIT_FAILURE);
  }
  if((ratio> min_ && ratio< max_) || (ratio>min_ && max_ ==-1) || (ratio<max_ && min_ ==-1)) return true;
  return false;

}

PTWhadCut::PTWhadCut(Context & ctx, float min, float max, const std::string & hyp_name):min_(min), max_(max){
 recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
}

bool PTWhadCut::passes(const uhh2::Event & event){
  BprimeContainer hyp = event.get(recohyp);
  if(hyp.get_wHad().pt()>min_ && (hyp.get_wHad().pt() <max_ || max_==-1))return true;
  return false;
}

ForwardJetPtEtaCut::ForwardJetPtEtaCut(float minEta, float maxEta, float minPt, float maxPt):minEta_(minEta),maxEta_(maxEta),minPt_(minPt),maxPt_(maxPt){}

bool ForwardJetPtEtaCut::passes(const uhh2::Event & event){
  vector<Jet> jets = *event.jets;
  Jet maxEtaJet = jets.at(0);
  for(auto jet : jets)
    if(fabs(jet.eta())> fabs(maxEtaJet.eta()) && minPt_<jet.pt() && maxPt_ > jet.pt()) maxEtaJet = jet;
  return fabs(maxEtaJet.eta())>minEta_ && maxEtaJet.pt()>minPt_ && (maxEtaJet.eta() < maxEta_ || maxEta_==-1) && (maxEtaJet.pt()<maxPt_ || maxPt_==-1);
}
