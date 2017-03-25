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

HTLepSelection::HTLepSelection(Context & ctx, double HTLepmin_, std::string metname_): HTLepmin(HTLepmin_),metname(metname_){
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  if(!metname.empty())
    h_met = ctx.get_handle<MET>(metname);
}

bool HTLepSelection::passes(const Event & event){
  if(metname.empty()){
    return HTLepmin < event.met->pt()+event.get(h_primlep).pt();
  }
  else{
    //cout<<"HTLep "<<event.get(h_met).pt()+event.get(h_primlep).pt() <<endl;
    return HTLepmin < event.get(h_met).pt()+event.get(h_primlep).pt();
  }
}
METSelection::METSelection(double METmin_): METmin(METmin_){
  metcoll="";
}
METSelection::METSelection(double METmin_, Context & ctx, string metcoll_): METmin(METmin_),metcoll(metcoll_){
  if(!metcoll.empty())
    h_met = ctx.get_handle<MET>(metcoll);
}

bool METSelection::passes(const Event & event){
  if(metcoll.empty())
    return METmin < event.met->pt();
  else{
    //cout<<"MET "<< event.get(h_met).pt() <<endl;
    return METmin < event.get(h_met).pt();
  }
}

bool RelIso::passes(const Event & event){
  assert((event.muons || event.electrons));
  Muon leading_muon;
  for(auto & muon : *event.muons)
    if(muon.pt()>leading_muon.pt()) leading_muon = muon;
  Electron leading_ele;
  for(auto & ele : *event.electrons)
    if(ele.pt()>leading_ele.pt())leading_ele = ele;
  if(event.muons->size()){
    return leading_muon.relIso() < reliso;
  }
  else
    return leading_ele.relIso() < reliso;
}
  
TwoDCut::TwoDCut(uhh2::Context & ctx, std::string jetcoll_,float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {
  h_jets = ctx.get_handle<std::vector<Jet> >(jetcoll_);
}

bool TwoDCut::passes(const Event & event){
  vector<Jet> all_jets;
  if(jetcoll.empty())
    all_jets = *event.jets;
  else 
    all_jets = event.get(h_jets);

  if(all_jets.size()==0)return false;
  //assert((event.muons || event.electrons));

  Muon leading_muon;
  for(auto & muon : *event.muons)
    if(muon.pt()>leading_muon.pt()) leading_muon = muon;

  Electron leading_ele;
  for(auto & ele : *event.electrons)
    if(ele.pt()>leading_ele.pt())leading_ele = ele;

  vector<Jet> jets;
  for(const auto jet :all_jets)
    if(abs(jet.eta())<2.4 && jet.pt()>=30) jets.push_back(jet);
  if(jets.size()==0)return false;

  float drmin, ptrel;  
  if(event.muons->size()){
    std::tie(drmin, ptrel) = drmin_pTrel(leading_muon, jets);
  }
  else{
    std::tie(drmin, ptrel) = drmin_pTrel(leading_ele, jets);
  }
  /*
  Jet close_jet = jets[0];
  for(auto jet : jets){
    if(deltaR(leading_muon.v4(),close_jet.v4()) > deltaR(leading_muon.v4(),jet.v4()))
       close_jet = jet;
  }
  cout<<"======================"<<endl;
  cout<<"Jet"<<endl;
  cout<<"pt "<<close_jet.pt()<<" eta "<<close_jet.eta()<<endl;
  cout<<"Muon"<<endl;
  cout<<"pt "<<leading_muon.pt()<<" eta "<<leading_muon.eta()<<endl;
  cout<<"ptrel "<<ptrel<<" drmin "<<drmin<<endl;
  */

  //cout<<"drmin "<<drmin<<" ptrel "<<ptrel <<endl;

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

ForwardJetPtEtaCut::ForwardJetPtEtaCut(Context & ctx,float minEta, float maxEta, float minPt, float maxPt, float minDrmin, float minEnergy, std::string hyp_name):minEta_(minEta),maxEta_(maxEta),minPt_(minPt),maxPt_(maxPt),minDrmin_(minDrmin),minEnergy_(minEnergy),hyp_name_(hyp_name){
  if(!hyp_name.empty())
    recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
}

bool ForwardJetPtEtaCut::passes(const uhh2::Event & event){
  vector<Jet> jets = *event.jets;
  LorentzVector maxEtaJet(0,0,0,0);
  float minimaldr = -1;
  if(hyp_name_.empty()){
    for(auto jet : jets){
      if(fabs(jet.eta())> fabs(maxEtaJet.eta()) && minPt_<jet.pt() && (maxPt_ > jet.pt() || maxPt_ ==-1) ) 
	maxEtaJet = jet.v4();
    }
  }
  else{
    BprimeContainer hyp = event.get(recohyp);
    maxEtaJet=hyp.get_forwardJet();
  }
  for(auto jet : jets){
    if((deltaR(maxEtaJet,jet.v4()) < minimaldr || minimaldr==-1 ) && deltaR(maxEtaJet,jet.v4()) > 0) 
      minimaldr = deltaR(maxEtaJet,jet.v4());
  }
   


  //cout<<"E "<<maxEtaJet.v4().E()<<" > "<<minEnergy_<<" min dR "<<minimaldr<<" > "<<minDrmin_<<endl;
  return fabs(maxEtaJet.eta())>= minEta_ && (fabs(maxEtaJet.eta()) <= maxEta_ || maxEta_==-1) && maxEtaJet.energy() >= minEnergy_;
  //return fabs(maxEtaJet.eta())> minEta_ && (fabs(maxEtaJet.eta()) < maxEta_ || maxEta_==-1) && maxEtaJet.energy() > minEnergy_ && minimaldr > minDrmin_;
  //return ((fabs(maxEtaJet.eta())>= minEta_  && (fabs(maxEtaJet.eta()) <= maxEta_ || maxEta_==-1)) || (minimaldr >= minDrmin_ && minDrmin_ >=0) ) && (maxEtaJet.energy() >= minEnergy_ || minEnergy_ ==0.);
  
}


NSubJetCut::NSubJetCut(int min_subjets_, int max_subjets_, int first_topjet_ , int last_topjet_):min_subjets(min_subjets_),max_subjets(max_subjets_),first_topjet(first_topjet_),last_topjet(last_topjet_){}

bool NSubJetCut::passes(const uhh2::Event & event){
  vector<TopJet> topjets = *event.topjets;
  bool result = false;
  for(unsigned int i =0; i< topjets.size(); i++){
    if(i+1>= first_topjet && (i+1<= last_topjet || last_topjet ==-1)){
      if(topjets.at(i).subjets().size() >=  min_subjets && (topjets.at(i).subjets().size()<= max_subjets || max_subjets ==-1 ))
	{
	result=true;
      }
      else
	return false;
    } 
  }
  return result;
}


TopJetMassCut::TopJetMassCut(float mass_): mass(mass_){}

bool TopJetMassCut::passes(const uhh2::Event & event){
  TopJet topjet = (*event.topjets)[0];
  return topjet.prunedmass() > mass;
}
