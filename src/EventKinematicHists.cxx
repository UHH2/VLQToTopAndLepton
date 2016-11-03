#include "UHH2/VLQToTopAndLepton/include/EventKinematicHists.h"

using namespace std;
using namespace uhh2;


EventKinematicHists::BaseHists EventKinematicHists::book_BaseHists(const std::string & name, const std::string & label, double minMass, double maxMass, double minPt, double maxPt, double minE, double maxE){
  BaseHists hists;
  hists.pt     = book<TH1F>("pt_"+name,"p_{T} "+ label+" [GeV]",100,minPt,maxPt);
  hists.eta    = book<TH1F>("eta_"+name,"#eta "+label,100,-5,5);
  hists.phi    = book<TH1F>("phi_"+name,"#phi "+label,100,-3.2,3.2);
  hists.mass   = book<TH1F>("mass_"+name,"Mass "+label+" [GeV]",30,minMass,maxMass);
  hists.energy = book<TH1F>("energy_"+name,"E "+label+" [GeV]",30,minE,maxE);
  return hists;
}
template<typename T>
void EventKinematicHists::fill_BaseHists(const T & particle, BaseHists & hists, double weight){
  hists.pt->Fill(particle.pt(),weight);
  hists.eta->Fill(particle.eta(),weight);
  hists.phi->Fill(particle.phi(),weight);
  hists.mass->Fill(sqrt(particle.M2()),weight);
  hists.energy->Fill(particle.E(),weight);
}

EventKinematicHists::EventKinematicHists(uhh2::Context & ctx, const std::string & dirname, std::string hyp_name): Hists(ctx, dirname){
  primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  leadingAk4_lepton_dr = book<TH1F>("leadingAk4_lepton_dr","#Delta R (AK4,lepton)", 100, 0, 8);
  leadingAk4_lepton_dphi  = book<TH1F>("leadingAk4_lepton_dphi","#Delta #Phi (Ak4,lepton)", 100, 0, 3.41);
  leadingAk8_lepton_dr = book<TH1F>("leadingAk8_lepton_dr","#Delta R (Ak8,lepton)", 100, 0, 8);
  leadingAk8_lepton_dphi  = book<TH1F>("leadingAk8_lepton_dphi","#Delta #Phi (Ak8,lepton)", 100, 0, 3.41); 
  massleadingAk8_lepton_dr = book<TH1F>("massleadingAk8_lepton_dr","mass leading Ak8 #Delta R (Ak8,lepton)", 100, 0, 8);
  massleadingAk8_lepton_dphi  = book<TH1F>("massleadingAk8_lepton_dphi","mass leading Ak8 #Delta #Phi (Ak8,lepton)", 100, 0, 3.41); 
  massleadingAk8_prunedmass  = book<TH1F>("massleadingAk8_prunedmass","mass leading Ak8 pruned mas", 100, 0, 500); 
  leadingAk4_lepton_dphi_mass  = book<TH2F>("leadingAk4_lepton_dphi_mass","#Delta #Phi(Ak4,lep) mass Ak4", 100, 0, 3.41, 100, 0, 500); 
  leadingEtadrmin_eta = book<TH2F>("leadingEtadrmin_eta","leadingEtadrmin #eta", 100, 0, 5, 100, -5, 5); 
  leadingEtadrmin_energy = book<TH2F>("leadingEtadrmin_energy","leadingEtadrmin Energy", 100, 0, 5, 100, 0, 2000); 
  energy_eta = book<TH2F>("energy_eta","Energy #eta", 100, 0, 2000, 100, -5, 5); 
  centrality_ak4 = book<TH1F>("centrality_ak4","Centrality Ak4", 100, 0, 1); 
  centrality_ak8 = book<TH1F>("centrality_ak8","Centrality Ak8", 100, 0, 1); 
  leadingAk8Hists = book_BaseHists("leadingAk8Jet", "leading p_{T} Ak8", 0, 400,0, 1200);
  leadingAk4Hists = book_BaseHists("leadingAk4Jet", "leading p_{T} Ak4", 0, 400,0, 1200);
  leadingEtaJetHists = book_BaseHists("leadingEtaJet", "leading #eta -Jet", 0, 300,0, 1000);
  leadingEtadrmin = book<TH1F>("leadingEtadrmin","#Delta R_{min} (leading eta Ak4,Ak4) leading eta jet", 100, 0, 5);
  secondEtaJetHists  = book_BaseHists("secondEtaJet", "second #eta -Jet", 0, 300,0, 1000);
  recoForward = book_BaseHists("RecoForward", "reco forward-Jet", 0, 300,0, 1000);
  forwardjet_iso = book<TH1F>("iso_recoforward","iso for the reconstructed forward jet", 200, 0, 10);
  noforward_Jet = book<TH1F>("number_noforward","No Forward jet in the Event",5,0.5,1.5);
  lbHists = book_BaseHists("mass_lb","lep+b [GeV]",0,2000.);
  METlHists = book_BaseHists("METl","MET+lep [GeV]",0,2000.);
  hyp_name_= hyp_name;
  ///cout<<hyp_name_<<" "<< hyp_name<<endl;
  if(!hyp_name_.empty())
    recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
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

  if(!hyp_name_.empty()){
    BprimeContainer hyp = event.get(recohyp);
    //cout<<"pt "<<hyp.get_forwardJet().pt()<<" phi "<<hyp.get_forwardJet().phi()<<" eta "<<hyp.get_forwardJet().eta()<<" E "<<hyp.get_forwardJet().E()<<endl;
    if(hyp.get_forwardJet().pt()>0){
      fill_BaseHists(hyp.get_forwardJet(),recoForward,weight);
      double iso_new = 999999999;
      for(auto & ak4 : *event.jets){
	if(iso_new > deltaR(ak4.v4(),hyp.get_forwardJet()) && deltaR(ak4.v4(),hyp.get_forwardJet())>0)
	  iso_new = deltaR(ak4.v4(),hyp.get_forwardJet());
      }
      forwardjet_iso->Fill(iso_new,weight);
    }
    else{
      noforward_Jet->Fill(1,weight);
    }
  }

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

    double drmin_leadingeta_ak4 = -1;
    for(auto & ak4 : *event.jets){
      if((deltaR(ak4,leadingEtaJet) < drmin_leadingeta_ak4 || drmin_leadingeta_ak4 ==-1) && deltaR(ak4,leadingEtaJet) > 0)
	drmin_leadingeta_ak4 = deltaR(ak4,leadingEtaJet);
    }
    leadingEtadrmin->Fill(drmin_leadingeta_ak4,weight);
    leadingEtadrmin_eta->Fill(drmin_leadingeta_ak4,leadingEtaJet.eta(),weight);
    leadingEtadrmin_energy->Fill(drmin_leadingeta_ak4,leadingEtaJet.E(),weight);
    energy_eta->Fill(leadingEtaJet.E(),leadingEtaJet.eta(),weight);
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


