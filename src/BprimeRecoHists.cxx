#include "UHH2/VLQToTopAndLepton/include/BprimeRecoHists.h"

using namespace std;
using namespace uhh2;

BprimeRecoHists::BaseHists BprimeRecoHists::book_BaseHists(const std::string & name, const std::string & label, double minPt, double maxPt){
  BaseHists hists;
  hists.pt   = book<TH1F>("pt_"+name,"p_{T} "+label,100,minPt,maxPt);
  hists.eta  = book<TH1F>("eta_"+name,"#eta "+label,100,-4,4);
  hists.phi  = book<TH1F>("phi_"+name,"#phi "+label,100,-3.2,3.2);
  hists.mass = book<TH1F>("phi_"+name,"#phi "+label,100,0,600);
}

template<typename T>
void BprimeRecoHists::fill_BaseHists(const T & particle, BaseHists & hists, double weight){
  hists.pt->Fill(particle.pt(),weight);
  hists.eta->Fill(particle.eta(),weight);
  hists.phi->Fill(particle.phi(),weight);
  hists.mass->Fill(sqrt(particle.M2()),weight);
}

BprimeRecoHists::BprimeRecoHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  wHad_all = book_BaseHists("wHad_all","all W_{had} Hypothesis");
  wLep_all = book_BaseHists("wLep_all","all W_{lep} Hypothesis"); 
  topLep_all = book_BaseHists("topLep_all","all Top_{lep} Hypothesis");
  topHad_all = book_BaseHists("topHad_all","all Top_{had} Hypothesis");
  
  deltaR_w_all    = book<TH1F>("deltaR_w","#Delta R (W_{lep},W_{had})", 100, 0, 8);
  deltaPhi_w_all  = book<TH1F>("deltaPhi_w","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);
  
  wHad_best = book_BaseHists("wHad_best","best W_{had} Hypothesis");
  wLep_best = book_BaseHists("wLep_best","best W_{lep} Hypothesis"); 
  topLep_best = book_BaseHists("topLep_best","best Top_{lep} Hypothesis");
  topHad_best = book_BaseHists("topHad_best","best Top_{had} Hypothesis");

  deltaR_w_best    = book<TH1F>("deltaR_w_best","#Delta R (W_{lep},W_{had})", 100, 0, 8);
  deltaPhi_w_best  = book<TH1F>("deltaPhi_w_best","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);

  hyps = ctx.get_handle<std::vector<BprimeContainer>>("BprimeReco");

}

BprimeRecoHists::~BprimeRecoHists(){}


void BprimeRecoHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  
  LorentzVector whad_best(0,0,0,0);
  LorentzVector wlep_best(0,0,0,0);
  LorentzVector topJets_best(0,0,0,0);

  for(auto hyp :  event.get(hyps)){
    LorentzVector whad = hyp.get_wHad();
    LorentzVector wlep = hyp.get_wLep();
    LorentzVector topJets = hyp.get_topJets();

    fill_BaseHists(whad, wHad_all, weight);
    fill_BaseHists(wlep, wLep_all, weight);
    fill_BaseHists(topJets+whad,topHad_all,weight);
    fill_BaseHists(topJets+wlep,topLep_all,weight);

    deltaR_w_all->Fill(deltaR(whad,wlep),weight);
    deltaPhi_w_all->Fill(deltaPhi(whad,wlep),weight);

    if(abs(whad_best.M()-82) > abs(whad.M()-82)){
      whad_best=whad;
      wlep_best=wlep;
      topJets_best=topJets;
    }
  }

  whad_best_pt->Fill(whad_best.pt(), weight);
  fill_BaseHists(whad_best, wHad_best, weight);
  fill_BaseHists(wlep_best, wLep_best, weight);
  fill_BaseHists(topJets_best+whad_best,topHad_best,weight);
  fill_BaseHists(topJets_best+wlep_best,topLep_best,weight);

  deltaR_w_best->Fill(deltaR(whad_best,wlep_best),weight);
  deltaPhi_w_best->Fill(deltaPhi(whad_best,wlep_best),weight);
}
