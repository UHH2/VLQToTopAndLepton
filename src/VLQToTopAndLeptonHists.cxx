#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

VLQToTopAndLeptonHists::VLQToTopAndLeptonHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  N_jets = book<TH1F>("N_jets", "N_{jets}", 20, 0, 20); 
  eta_jet1 = book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -5, 5);
  eta_jet2 = book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -5, 5);
  eta_jet3 = book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -5, 5);
  eta_jet4 = book<TH1F>("eta_jet4", "#eta^{jet 4}", 40, -5, 5);

  pt_eta_jet1 = book<TH2F>("pt_eta_jet1", "x: p_{T}, y: #eta jet 1", 100, 0, 150, 100, 0, 5);
  pt_eta_jet2 = book<TH2F>("pt_eta_jet2", "x: p_{T}, y: #eta jet 2", 100, 0, 150, 100, 0, 5);
  pt_eta_jet3 = book<TH2F>("pt_eta_jet3", "x: p_{T}, y: #eta jet 3", 100, 0, 150, 100, 0, 5);
  pt_eta_jet4 = book<TH2F>("pt_eta_jet4", "x: p_{T}, y: #eta jet 4", 100, 0, 150, 100, 0, 5);
  
  // leptons
  N_mu = book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  pt_mu = book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  eta_mu = book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.4, 2.4);
  reliso_mu = book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  // primary vertices
  N_pv = book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);
}


void VLQToTopAndLeptonHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  N_jets->Fill(Njets, weight);
  
  if(Njets>=1){
    eta_jet1->Fill(jets->at(0).eta(), weight);
    pt_eta_jet1->Fill(jets->at(0).pt(),fabs(jets->at(0).eta()),weight);
  }
  if(Njets>=2){
    eta_jet2->Fill(jets->at(1).eta(), weight);
    pt_eta_jet2->Fill(jets->at(1).pt(),fabs(jets->at(1).eta()),weight);
  }
  if(Njets>=3){
    eta_jet3->Fill(jets->at(2).eta(), weight);
    pt_eta_jet3->Fill(jets->at(2).pt(),fabs(jets->at(2).eta()),weight);
  }
  if(Njets>=4){
    eta_jet4->Fill(jets->at(3).eta(), weight);
    pt_eta_jet4->Fill(jets->at(3).pt(),fabs(jets->at(3).eta()),weight);
  }

  

  int Nmuons = event.muons->size();
  N_mu->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      pt_mu->Fill(thismu.pt(), weight);
      eta_mu->Fill(thismu.eta(), weight);
      reliso_mu->Fill(thismu.relIso(), weight);
  }
  
  int Npvs = event.pvs->size();
  N_pv->Fill(Npvs, weight);
}

VLQToTopAndLeptonHists::~VLQToTopAndLeptonHists(){}
