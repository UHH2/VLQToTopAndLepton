#include "UHH2/VLQToTopAndLepton/include/BprimeRecoHists.h"

using namespace std;
using namespace uhh2;

BprimeRecoHists::BaseHists BprimeRecoHists::book_BaseHists(const std::string & name, const std::string & label, double minMass, double maxMass, double minPt, double maxPt){
  BaseHists hists;
  hists.pt   = book<TH1F>("pt_"+name,"p_{T} "+label,100,minPt,maxPt);
  hists.eta  = book<TH1F>("eta_"+name,"#eta "+label,100,-4,4);
  hists.phi  = book<TH1F>("phi_"+name,"#phi "+label,100,-3.2,3.2);
  hists.mass = book<TH1F>("mass_"+name,"Mass "+label,50,minMass,maxMass);
  return hists;
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
  Bprime = book_BaseHists("Bprime","Bprime Hypothesis",50,3000); 
  Bprime_lep = book_BaseHists("Bprime_lep","Bprime lep. Hypothesis",50,3000); 
  Bprime_had = book_BaseHists("Bprime_had","Bprime had. Hypothesis",50,3000); 

  deltaR_w_all    = book<TH1F>("deltaR_w","#Delta R (W_{lep},W_{had})", 100, 0, 8);
  deltaPhi_w_all  = book<TH1F>("deltaPhi_w","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);
  
  wHad_best = book_BaseHists("wHad_best","best W_{had} Hypothesis");
  wLep_best = book_BaseHists("wLep_best","best W_{lep} Hypothesis"); 
  topLep_best = book_BaseHists("topLep_best","best Top_{lep} Hypothesis");
  topHad_best = book_BaseHists("topHad_best","best Top_{had} Hypothesis");

  deltaR_w_best    = book<TH1F>("deltaR_w_best","#Delta R (W_{lep},W_{had})", 100, 0, 4);
  deltaPhi_w_best  = book<TH1F>("deltaPhi_w_best","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);

  recotype_h = book<TH1F>("recoType","0 hadronic Top, 1 leptonic Top",2, -0.5, 1.5);
  chiDis = book<TH1F>("chiVal","#Chi^{2} Value",100, 0,60);
  chiDis_lep = book<TH1F>("chiVal_lep","#Chi^{2} Value lep. Top",100, 0,60);
  chiDis_had = book<TH1F>("chiVal_had","#Chi^{2} Value had. Top",100, 0,60);
  ttbar_chi =  book<TH1F>("ttbarChi","#Chi^{2} Value for ttbar",100, 0,100);

  wHad_res_pt   = book<TH1F>("wHad_res_pt","W_{had} Resolution p_{T}", 100, -3, 2);
  wHad_res_E    = book<TH1F>("wHad_res_E","W_{had} Resolution E", 100, -3, 2); 
  wHad_res_mass = book<TH1F>("wHad_res_mass","W_{had} Resolution Mass", 100, -3, 2); 
  wHad_res_phi  = book<TH1F>("wHad_res_phi","W_{had} Resolution #phi", 100, 0, 3.14); 
  wHad_res_eta  = book<TH1F>(" wHad_res_eta","W_{had} Resolution #eta", 100, 0, 4); 
  wHad_res_deltaR  = book<TH1F>(" wHad_res_deltaR","W_{had} Resolution #Delta R", 100, 0, 4); 
  wLep_res_pt   = book<TH1F>("wLep_res_pt","W_{lep} Resolution p_{T}", 100, -3, 2);
  wLep_res_E    = book<TH1F>("wLep_res_E","W_{lep} Resolution E", 100, -3, 2); 
  wLep_res_mass = book<TH1F>("wLep_res_mass","W_{lep} Resolution Mass", 100, -3, 2); 
  wLep_res_phi  = book<TH1F>("wLep_res_phi","W_{lep} Resolution #phi", 100, 0, 3.14); 
  wLep_res_eta  = book<TH1F>("wLep_res_eta","W_{lep} Resolution #eta", 100, 0., 4); 
  wLep_res_deltaR  = book<TH1F>("wLep_res_deltaR","W_{lep} Resolution #Delta R", 100, 0., 4); 

  topHad_res_pt   = book<TH1F>("topHad_res_pt","Top_{had} Resolution p_{T}", 100, -3, 2);
  topHad_res_E    = book<TH1F>("topHad_res_E","Top_{had} Resolution E", 100, -3, 2); 
  topHad_res_mass = book<TH1F>("topHad_res_mass","Top_{had} Resolution Mass", 100, -3, 2); 
  topHad_res_phi  = book<TH1F>("topHad_res_phi","Top_{had} Resolution #phi", 100, 0, 3.14); 
  topHad_res_eta  = book<TH1F>("topHad_res_eta","Top_{had} Resolution #eta", 100, 0, 4); 
  topHad_res_deltaR  = book<TH1F>("topHad_res_deltaR","Top_{had} Resolution #Delta R", 100, 0, 4); 
  topLep_res_pt   = book<TH1F>("topLep_res_pt","Top_{lep} Resolution p_{T}", 100, -3, 2);
  topLep_res_E    = book<TH1F>("topLep_res_E","Top_{lep} Resolution E", 100, -3, 2); 
  topLep_res_mass = book<TH1F>("topLep_res_mass","Top_{lep} Resolution Mass", 100, -3, 2); 
  topLep_res_phi  = book<TH1F>("topLep_res_phi","Top_{lep} Resolution #phi", 100, 0, 3.14); 
  topLep_res_eta  = book<TH1F>("topLep_res_eta","Top_{lep} Resolution #eta", 100, 0., 4); 
  topLep_res_deltaR  = book<TH1F>("topLep_res_deltaR","Top_{lep} Resolution #Delta R", 100, 0., 4); 

  topReco_dR_pT_lep    = book<TH2F>("topReco_dR_pT_lep","Gen and Reco Top #Delta R pT", 100, 0, 2, 100,0,800);  
  topReco_dR_pTres_lep = book<TH2F>("topReco_dR_pTres_lep","Gen and Reco Top #Delta R pT", 100, 0, 2, 100,-2,2);  
  topReco_dR_pT_had    = book<TH2F>("topReco_dR_pT_had","Gen and Reco Top #Delta R pT", 100, 0, 2, 100,0,800);  
  topReco_dR_pTres_had = book<TH2F>("topReco_dR_pTres_had","Gen and Reco Top #Delta R pT", 100, 0, 2, 100,-2,2);  
  wReco_dR_pTres_lep   = book<TH2F>("wReco_dR_pTres_lep","Gen and Reco W #Delta R / pT", 100, 0, 2, 100,-2,2);  
  wReco_dR_pT_lep      = book<TH2F>("wReco_dR_pT_lep","Gen and Reco W #Delta R / pT", 100, 0, 2, 100,0,800);  
  wReco_dR_pTres_had   = book<TH2F>("wReco_dR_pTres_had","Gen and Reco W #delta R / pT", 100, 0, 2, 100,-2,2);  
  wReco_dR_pT_had      = book<TH2F>("wReco_dR_pT_had","Gen and Reco W #Delta R / pT", 100, 0, 2, 100,0,800);  

  hyps = ctx.get_handle<std::vector<BprimeContainer>>("BprimeReco");
  gen = ctx.get_handle<BprimeGenContainer>("BprimeGen");
}

BprimeRecoHists::~BprimeRecoHists(){}


void BprimeRecoHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  
  LorentzVector whad_best(0,0,0,0);
  LorentzVector wlep_best(0,0,0,0);
  LorentzVector topJets_best(0,0,0,0);
  BprimeGenContainer GenInfo = event.get(gen);
  double chiVal = 99999999999.;
  double ttbarchi = 99999999999.;
  int recotype = 0;

  for(auto hyp :  event.get(hyps)){
    LorentzVector whad = hyp.get_wHad();
    LorentzVector wlep = hyp.get_wLep();
    LorentzVector topJets = hyp.get_topJets();

    fill_BaseHists(whad, wHad_all, weight);
    fill_BaseHists(wlep, wLep_all, weight);

    if(fabs((whad).M()-175)<fabs((topJets+wlep).M()-175)){
      fill_BaseHists(topJets+whad,topHad_all,weight);
    }
    else
      fill_BaseHists(topJets+wlep,topLep_all,weight);

    deltaR_w_all->Fill(deltaR(whad,wlep),weight);
    deltaPhi_w_all->Fill(deltaPhi(whad,wlep),weight);

    double ttbar = ((topJets+wlep).M()-174)*((topJets+wlep).M()-174)/(18*18)+(whad.M()-181)*(whad.M()-181)/(15*15);
    if(ttbar < ttbarchi)ttbarchi=ttbar;
  
    double tophad_chi = ((topJets+whad).M()-178)*((topJets+whad).M()-178)/(14*14)+(whad.M()-82)*(whad.M()-82)/(20*20);
    double toplep_chi = ((topJets+wlep).M()-175)*((topJets+wlep).M()-175)/(14*14)+(whad.M()-82)*(whad.M()-82)/(20*20);
    if(chiVal > tophad_chi || chiVal> toplep_chi){
    //if(fabs((topJets_best+whad_best+wlep_best).M()-1200) > fabs((topJets+whad+wlep).M()-1200)){
      whad_best=whad;
      wlep_best=wlep;
      topJets_best=topJets;
      chiVal = tophad_chi > toplep_chi ? toplep_chi : tophad_chi;
      recotype = tophad_chi > toplep_chi ? 0 : 1;
    }
  }
  
  if(recotype==0){
    if((topJets_best+wlep_best).M()<100) return;
  }
  //if(ttbarchi<10)return; 
  fill_BaseHists(whad_best, wHad_best, weight);
  fill_BaseHists(wlep_best, wLep_best, weight);
  if(recotype==1){
    fill_BaseHists(topJets_best+whad_best,topHad_best,weight);
    fill_BaseHists(topJets_best+whad_best+wlep_best,Bprime_had,weight);
    chiDis_lep->Fill(chiVal,weight);
  }
  else{
    fill_BaseHists(topJets_best+wlep_best,topLep_best,weight);
    fill_BaseHists(topJets_best+whad_best+wlep_best,Bprime_lep,weight);
    chiDis_had->Fill(chiVal,weight);
  }

  

  if(GenInfo.get_wHad().pt()>0){
    wHad_res_pt->Fill((GenInfo.get_wHad().pt()-whad_best.pt())/GenInfo.get_wHad().pt(),weight);   
    wHad_res_E->Fill((GenInfo.get_wHad().E()-whad_best.E())/GenInfo.get_wHad().E(),weight);       
    wHad_res_mass->Fill((GenInfo.get_wHad().M()-whad_best.M())/GenInfo.get_wHad().M(),weight);    
    wHad_res_phi->Fill(deltaPhi(GenInfo.get_wHad(),whad_best),weight);     
    wHad_res_eta->Fill(abs(GenInfo.get_wHad().eta()-whad_best.eta()),weight);     
    wHad_res_deltaR->Fill(deltaR(GenInfo.get_wHad(),whad_best),weight);
    wReco_dR_pTres_had->Fill(deltaR(GenInfo.get_wHad(),whad_best),(GenInfo.get_wHad().pt()-whad_best.pt())/GenInfo.get_wHad().pt(),weight);
    wReco_dR_pT_had->Fill(deltaR(GenInfo.get_wHad(),whad_best),GenInfo.get_wHad().pt(),weight);     
  }  
  if(GenInfo.get_wLep().pt()>0){
    wLep_res_pt->Fill((GenInfo.get_wLep().pt()-wlep_best.pt())/GenInfo.get_wLep().pt(),weight);     
    wLep_res_E->Fill((GenInfo.get_wLep().E()-wlep_best.E())/GenInfo.get_wLep().E(),weight);       
    wLep_res_mass->Fill((GenInfo.get_wLep().M()-wlep_best.M())/GenInfo.get_wLep().M(),weight);    
    wLep_res_phi->Fill(deltaPhi(GenInfo.get_wLep(),wlep_best),weight);     
    wLep_res_eta->Fill(abs(GenInfo.get_wLep().eta()-wlep_best.eta()),weight);     
    wLep_res_deltaR->Fill(deltaR(GenInfo.get_wLep(),wlep_best),weight);
    wReco_dR_pTres_lep->Fill(deltaR(GenInfo.get_wLep(),wlep_best),(GenInfo.get_wLep().pt()-wlep_best.pt())/GenInfo.get_wLep().pt(),weight);  
    wReco_dR_pT_lep->Fill(deltaR(GenInfo.get_wLep(),wlep_best),GenInfo.get_wLep().pt(),weight);       
  }



  if(GenInfo.get_topLep().pt()>0 && recotype==0){
    topReco_dR_pT_lep->Fill(deltaR(GenInfo.get_topLep(),wlep_best+topJets_best),GenInfo.get_topLep().pt(),weight);   
    topReco_dR_pTres_lep->Fill(deltaR(GenInfo.get_topLep(),wlep_best+topJets_best),(GenInfo.get_topLep().pt()-(wlep_best+topJets_best).pt())/GenInfo.get_topLep().pt(),weight);
    topLep_res_pt->Fill((GenInfo.get_topLep().pt()-(wlep_best+topJets_best).pt())/GenInfo.get_topLep().pt(),weight);     
    topLep_res_E->Fill((GenInfo.get_topLep().E()-(wlep_best+topJets_best).E())/GenInfo.get_topLep().E(),weight);       
    topLep_res_mass->Fill((GenInfo.get_topLep().M()-(wlep_best+topJets_best).M())/GenInfo.get_topLep().M(),weight);    
    topLep_res_phi->Fill(deltaPhi(GenInfo.get_topLep(),wlep_best+topJets_best),weight);     
    topLep_res_eta->Fill(abs(GenInfo.get_topLep().eta()-(wlep_best+topJets_best).eta()),weight);     
    topLep_res_deltaR->Fill(deltaR(GenInfo.get_topLep(),wlep_best+topJets_best),weight);
  }
  if(GenInfo.get_topHad().pt()>0 && recotype==1){
    topReco_dR_pT_had->Fill(deltaR(GenInfo.get_topHad(),whad_best+topJets_best),GenInfo.get_topHad().pt(),weight);   
    topReco_dR_pTres_had->Fill(deltaR(GenInfo.get_topHad(),whad_best+topJets_best),(GenInfo.get_topHad().pt()-(whad_best+topJets_best).pt())/GenInfo.get_topHad().pt(),weight);
    topHad_res_pt->Fill((GenInfo.get_topHad().pt()-(whad_best+topJets_best).pt())/GenInfo.get_topHad().pt(),weight);   
    topHad_res_E->Fill((GenInfo.get_topHad().E()-(whad_best+topJets_best).E())/GenInfo.get_topHad().E(),weight);       
    topHad_res_mass->Fill((GenInfo.get_topHad().M()-(whad_best+topJets_best).M())/GenInfo.get_topHad().M(),weight);    
    topHad_res_phi->Fill(deltaPhi(GenInfo.get_topHad(),whad_best+topJets_best),weight);     
    topHad_res_eta->Fill(abs(GenInfo.get_topHad().eta()-(whad_best+topJets_best).eta()),weight);     
    topHad_res_deltaR->Fill(deltaR(GenInfo.get_topHad(),whad_best+topJets_best),weight);
  }

  fill_BaseHists(topJets_best+whad_best+wlep_best,Bprime,weight);

  ttbar_chi->Fill(ttbarchi,weight);
  deltaR_w_best->Fill(deltaR(whad_best,wlep_best),weight);
  deltaPhi_w_best->Fill(deltaPhi(whad_best,wlep_best),weight);
  recotype_h->Fill(recotype,weight);
  chiDis->Fill(chiVal,weight);
}


//double calc_Chi(LorentzVector top, LorentzVector w){
//}
