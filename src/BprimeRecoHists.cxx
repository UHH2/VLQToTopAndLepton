#include "UHH2/VLQToTopAndLepton/include/BprimeRecoHists.h"

using namespace std;
using namespace uhh2;

BprimeRecoHists::BprimeRecoHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  whad_pt   =  book<TH1F>("whad_pt"  ,"had W p_T" , 100, 0, 1400);
  whad_phi  =  book<TH1F>("whad_phi" ,"had W #phi", 100, -3.2, 3.2);
  whad_eta  =  book<TH1F>("whad_eta" ,"had W #eta", 100, -4, 4);
  whad_mass =  book<TH1F>("whad_mass","had W Mass", 100, 0, 200);

  wlep_pt   =  book<TH1F>("wlep_pt"  ,"lep W p_T" , 100, 0, 1400);
  wlep_phi  =  book<TH1F>("wlep_phi" ,"lep W #phi", 100, -3.2, 3.2);
  wlep_eta  =  book<TH1F>("wlep_eta" ,"lep W #eta", 100, -4, 4);
  wlep_mass =  book<TH1F>("wlep_mass","lep W Mass", 100, 0, 250);
  
  deltaR_w    = book<TH1F>("deltaR_w","#Delta R (W_{lep},W_{had})", 100, 0, 8);
  deltaPhi_w  = book<TH1F>("deltaPhi_w","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);


  whad_best_pt   =  book<TH1F>("whad_best_pt"  ,"had best W p_T" , 100, 0, 1400);
  whad_best_phi  =  book<TH1F>("whad_best_phi" ,"had best W #phi", 100, -3.2, 3.2);
  whad_best_eta  =  book<TH1F>("whad_best_eta" ,"had best W #eta", 100, -4, 4);
  whad_best_mass =  book<TH1F>("whad_best_mass","had best W Mass", 100, 0, 200);

  wlep_best_pt   =  book<TH1F>("wlep_best_pt"  ,"lep best W p_T" , 100, 0, 1400);
  wlep_best_phi  =  book<TH1F>("wlep_best_phi" ,"lep best W #phi", 100, -3.2, 3.2);
  wlep_best_eta  =  book<TH1F>("wlep_best_eta" ,"lep best W #eta", 100, -4, 4);
  wlep_best_mass =  book<TH1F>("wlep_best_mass","lep best W Mass", 100, 0, 250);
  
  deltaR_w_best    = book<TH1F>("deltaR_w_best","#Delta R (W_{lep},W_{had})", 100, 0, 8);
  deltaPhi_w_best  = book<TH1F>("deltaPhi_w_best","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);




  hyps = ctx.get_handle<std::vector<BprimeContainer>>("BprimeReco");

}

BprimeRecoHists::~BprimeRecoHists(){}


void BprimeRecoHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  
  LorentzVector whad_best(0,0,0,0);
  LorentzVector wlep_best(0,0,0,0);


  for(auto hyp :  event.get(hyps)){
    LorentzVector whad = hyp.get_wHad();
    LorentzVector wlep = hyp.get_wLep();
    whad_pt->Fill(whad.pt(), weight);
    whad_phi->Fill(whad.phi(), weight);
    whad_eta->Fill(whad.eta(), weight);
    whad_mass->Fill(whad.M(), weight);
    wlep_pt->Fill(wlep.pt(), weight);
    wlep_phi->Fill(wlep.phi(), weight);
    wlep_eta->Fill(wlep.eta(), weight);
    wlep_mass->Fill(wlep.M(), weight);

    deltaR_w->Fill(deltaR(whad,wlep),weight);
    deltaPhi_w->Fill(deltaPhi(whad,wlep),weight);

    if(abs(whad_best.M()-82) > abs(whad.M()-82)){
      whad_best=whad;
      wlep_best=wlep;
    }
  }

  whad_best_pt->Fill(whad_best.pt(), weight);
  whad_best_phi->Fill(whad_best.phi(), weight);
  whad_best_eta->Fill(whad_best.eta(), weight);
  whad_best_mass->Fill(whad_best.M(), weight);
  wlep_best_pt->Fill(wlep_best.pt(), weight);
  wlep_best_phi->Fill(wlep_best.phi(), weight);
  wlep_best_eta->Fill(wlep_best.eta(), weight);
  wlep_best_mass->Fill(wlep_best.M(), weight);

  deltaR_w_best->Fill(deltaR(whad_best,wlep_best),weight);
  deltaPhi_w_best->Fill(deltaPhi(whad_best,wlep_best),weight);
}
