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
  leadingAk4_lepton_dr = book<TH1F>("leadingAk4_lepton_dr","#Delta R (AK4,lepton)", 100, 0, 8);
  leadingAk4_lepton_dphi  = book<TH1F>("leadingAk4_lepton_dphi","#Delta #Phi (Ak4,lepton)", 100, 0, 3.41);
  leadingAk8_lepton_dr = book<TH1F>("leadingAk8_lepton_dr","#Delta R (Ak8,lepton)", 100, 0, 8);
  leadingAk8_lepton_dphi  = book<TH1F>("leadingAk8_lepton_dphi","#Delta #Phi (Ak8,lepton)", 100, 0, 3.41); 
  massleadingAk8_lepton_dr = book<TH1F>("massleadingAk8_lepton_dr","mass leading Ak8 #Delta R (Ak8,lepton)", 100, 0, 8);
  massleadingAk8_lepton_dphi  = book<TH1F>("massleadingAk8_lepton_dphi","mass leading Ak8 #Delta #Phi (Ak8,lepton)", 100, 0, 3.41); 
  massleadingAk8_prunedmass  = book<TH1F>("massleadingAk8_prunedmass","mass leading Ak8 pruned mas", 100, 0, 500); 
  leadingAk4_lepton_dphi_mass  = book<TH2F>("leadingAk4_lepton_dphi_mass","#Delta #Phi(Ak4,lep) mass Ak4", 100, 0, 3.41, 100, 0, 500); 
  centrality_ak4 = book<TH1F>("centrality_ak4","Centrality Ak4", 100, 0, 1); 
  centrality_ak8 = book<TH1F>("centrality_ak8","Centrality Ak8", 100, 0, 1); 
  leadingAk8Hists = book_BaseHists("leadingAk8Jet", "leading p_{T} Ak8", 0, 400,0, 1200);
  leadingAk4Hists = book_BaseHists("leadingAk4Jet", "leading p_{T} Ak4", 0, 400,0, 1200);
  leadingEtaJetHists = book_BaseHists("leadingEtaJet", "leading #eta -Jet", 0, 300,0, 1000);
  secondEtaJetHists  = book_BaseHists("secondEtaJet", "second #eta -Jet", 0, 300,0, 1000);
  lbHists = book_BaseHists("mass_lb","lep+b [GeV]",0,2000.);
  METlHists = book_BaseHists("METl","MET+lep [GeV]",0,2000.);
}

EventKinematicHists::~EventKinematicHists(){}

void EventKinematicHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  LorentzVector lep = event.get(primlep).v4();
  fill_BaseHists(lep+event.met->v4(),METlHists,weight);
  double sum_ak4_pt = 0;
  double sum_ak4_E = 0;
  double sum_ak8_pt = 0;
  double sum_ak8_E =0;

  if(event.jets->size()>0 ){
    //cout<<event.jets->size()<<endl;
    LorentzVector leadingAk4 = (*event.jets)[0].v4();
    LorentzVector leadingEtaJet = (*event.jets)[0].v4();
    LorentzVector secondEtaJet(0,0,0,0);
    Jet highest_btag_score = (*event.jets)[0];
    for(auto & ak4 : *event.jets){
      sum_ak4_pt += ak4.pt();
      sum_ak4_E += ak4.v4().energy();
      if(ak4.pt()>leadingAk4.pt())
	leadingAk4 = ak4.v4();
      if(fabs(ak4.eta()) > fabs(leadingEtaJet.eta()))
	leadingEtaJet = ak4.v4();
      if(highest_btag_score.btag_combinedSecondaryVertex() < ak4.btag_combinedSecondaryVertex())
	highest_btag_score = ak4;
    }
    for(auto & ak4 : *event.jets)
      if(fabs(leadingEtaJet.eta())>fabs(ak4.eta()) && fabs(ak4.eta())>fabs(secondEtaJet.eta()))
	secondEtaJet= ak4.v4();
    leadingAk4_lepton_dr->Fill(deltaR(leadingAk4,lep),weight);
    leadingAk4_lepton_dphi->Fill(deltaPhi(leadingAk4,lep),weight);
    leadingAk4_lepton_dphi_mass->Fill(deltaPhi(leadingAk4,lep),leadingAk4.M(),weight);

    centrality_ak4->Fill(sum_ak4_pt/sum_ak4_E,weight);

    fill_BaseHists(leadingAk4, leadingAk4Hists, weight);
    fill_BaseHists(leadingEtaJet, leadingEtaJetHists, weight);
    fill_BaseHists(secondEtaJet, secondEtaJetHists, weight);
    fill_BaseHists(highest_btag_score.v4()+lep,lbHists,weight);
  }

  if(event.topjets->size()>0){
    LorentzVector leadingAk8 = (*event.topjets)[0].v4();
    TopJet massleadingAk8 = (*event.topjets)[0];
    for(auto & ak8 : *event.topjets){
      sum_ak8_pt += ak8.pt();
      sum_ak8_E += ak8.v4().energy();
      if(ak8.pt()>leadingAk8.pt())
	leadingAk8 = ak8.v4();
      if(ak8.prunedmass() > massleadingAk8.prunedmass())
	massleadingAk8 = ak8;
    }
    centrality_ak8->Fill(sum_ak8_pt/sum_ak8_E,weight);
    leadingAk8_lepton_dr->Fill(deltaR(leadingAk8,lep),weight);
    leadingAk8_lepton_dphi->Fill(deltaPhi(leadingAk8,lep),weight);
    massleadingAk8_lepton_dr->Fill(deltaR(massleadingAk8,lep),weight);
    massleadingAk8_lepton_dphi->Fill(deltaPhi(massleadingAk8,lep),weight);
    fill_BaseHists(leadingAk8, leadingAk8Hists, weight);
    massleadingAk8_prunedmass->Fill(massleadingAk8.prunedmass(),weight);
  }
}


