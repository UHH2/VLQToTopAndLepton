#include "UHH2/VLQToTopAndLepton/include/EventKinematicHists.h"

using namespace std;
using namespace uhh2;


EventKinematicHists::BaseHists EventKinematicHists::book_BaseHists(const std::string & name, const std::string & label, double minMass, double maxMass, double minPt, double maxPt){
  BaseHists hists;
  hists.pt   = book<TH1F>("pt_"+name,"p_{T} "+ label+" [GeV]",100,minPt,maxPt);
  hists.eta  = book<TH1F>("eta_"+name,"#eta "+label,100,-5,5);
  hists.phi  = book<TH1F>("phi_"+name,"#phi "+label,100,-3.2,3.2);
  hists.mass = book<TH1F>("mass_"+name,"Mass "+label+" [GeV]",30,minMass,maxMass);
  return hists;
}
template<typename T>
void EventKinematicHists::fill_BaseHists(const T & particle, BaseHists & hists, double weight){
  hists.pt->Fill(particle.pt(),weight);
  hists.eta->Fill(particle.eta(),weight);
  hists.phi->Fill(particle.phi(),weight);
  hists.mass->Fill(sqrt(particle.M2()),weight);
}

EventKinematicHists::EventKinematicHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  leadingAk8_lepton_dr = book<TH1F>("leadingAk8_lepton_dr","#Delta R (Ak8,lepton)", 100, 0, 8);
  leadingAk4_lepton_dr = book<TH1F>("leadingAk4_lepton_dr","#Delta R (AK4,lepton)", 100, 0, 8);
  leadingAk8_lepton_dphi  = book<TH1F>("leadingAk8_lepton_dphi","#Delta #Phi (Ak8,lepton)", 100, 0, 3.41); 
  leadingAk4_lepton_dphi  = book<TH1F>("leadingAk4_lepton_dphi","#Delta #Phi (Ak4,lepton)", 100, 0, 3.41);
  leadingAk8Hists = book_BaseHists("leadingAk8Jet", "leading p_{T} Ak8", 0, 400,0, 1200);
  leadingAk4Hists = book_BaseHists("leadingAk4Jet", "leading p_{T} Ak4", 0, 400,0, 1200);
  leadingEtaJetHists = book_BaseHists("leadingEtaJet", "leading #eta -Jet", 0, 300,0, 1000);
  secondEtaJetHists  = book_BaseHists("secondEtaJet", "second #eta -Jet", 0, 300,0, 1000);
}

EventKinematicHists::~EventKinematicHists(){}

void EventKinematicHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  LorentzVector lep = event.get(primlep).v4();
  if(event.jets->size()>0 ){
    //cout<<event.jets->size()<<endl;
    LorentzVector leadingAk4 = (*event.jets)[0].v4();
    LorentzVector leadingEtaJet = (*event.jets)[0].v4();
    LorentzVector secondEtaJet(0,0,0,0);
    for(auto & ak4 : *event.jets){
      if(ak4.pt()>leadingAk4.pt())
	leadingAk4 = ak4.v4();
      if(fabs(ak4.eta()) > fabs(leadingEtaJet.eta()))
	leadingEtaJet = ak4.v4();
    }
    for(auto & ak4 : *event.jets)
      if(fabs(leadingEtaJet.eta())>fabs(ak4.eta()) && fabs(ak4.eta())>fabs(secondEtaJet.eta()))
	secondEtaJet= ak4.v4();
    leadingAk4_lepton_dr->Fill(deltaR(leadingAk4,lep),weight);
    leadingAk4_lepton_dphi->Fill(deltaPhi(leadingAk4,lep),weight);
    fill_BaseHists(leadingAk4, leadingAk4Hists, weight);
    fill_BaseHists(leadingEtaJet, leadingEtaJetHists, weight);
    fill_BaseHists(secondEtaJet, secondEtaJetHists, weight);
  }

  if(event.topjets->size()>0){
    LorentzVector leadingAk8 = (*event.topjets)[0].v4();
    for(auto & ak8 : *event.topjets)
      if(ak8.pt()>leadingAk8.pt())
	leadingAk8 = ak8.v4();
    leadingAk8_lepton_dr->Fill(deltaR(leadingAk8,lep),weight);
    leadingAk8_lepton_dphi->Fill(deltaPhi(leadingAk8,lep),weight);
    fill_BaseHists(leadingAk8, leadingAk8Hists, weight);
  }
}


