//framework includes
#include "UHH2/VLQToTopAndLepton/include/Lepton2DHist.h"


//#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/core/include/Event.h"

//general includes
#include <iostream>

using namespace std;
using namespace uhh2;


Lepton2DHist::Lepton2DHist(Context & ctx, const string & dirname): Hists(ctx, dirname){
  leadingPT_muon_pt_eta      = book<TH2F>("leadingPT_muon_pt_eta"    , "leadingPT_muon_pt_eta", 300, 20, 800, 100, -2.4, 2.4);  
  leadingPT_electron_pt_eta  = book<TH2F>("leadingPT_electron_pt_eta", "leadingPT_muon_pt_eta", 300, 20, 800, 100, -2.4, 2.4);  
  secondPT_muon_pt_eta      = book<TH2F>("secondPT_muon_pt_eta"    , "secondPT_muon_pt_eta", 300, 20, 800, 100, -2.4, 2.4);  
  secondPT_electron_pt_eta  = book<TH2F>("secondPT_electron_pt_eta", "secondPT_muon_pt_eta", 300, 20, 800, 100, -2.4, 2.4);  
}


void Lepton2DHist::fill(const Event & event){
   double w = event.weight;
   
   if(event.muons->size()>0)
     leadingPT_muon_pt_eta->Fill(event.muons->at(0).pt(),event.muons->at(0).eta(),w);  
   if(event.electrons->size()>0)
     leadingPT_electron_pt_eta->Fill(event.electrons->at(0).pt(),event.electrons->at(0).eta(),w);  
   if(event.muons->size()>1)
     secondPT_muon_pt_eta->Fill(event.muons->at(1).pt(),event.muons->at(1).eta(),w);  
   if(event.electrons->size()>1)
     secondPT_electron_pt_eta->Fill(event.electrons->at(1).pt(),event.electrons->at(1).eta(),w);  


}


Lepton2DHist::~Lepton2DHist(){}


