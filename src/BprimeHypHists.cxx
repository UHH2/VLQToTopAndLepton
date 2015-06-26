#include "UHH2/VLQToTopAndLepton/include/BprimeHypHists.h"

using namespace std;
using namespace uhh2;

BprimeHypHists::BaseHists BprimeHypHists::book_BaseHists(const std::string & name, const std::string & label, double minMass, double maxMass, double minPt, double maxPt){
  BaseHists hists;
  hists.pt   = book<TH1F>("pt_"+name,"p_{T} "+label,100,minPt,maxPt);
  hists.eta  = book<TH1F>("eta_"+name,"#eta "+label,100,-4,4);
  hists.phi  = book<TH1F>("phi_"+name,"#phi "+label,100,-3.2,3.2);
  hists.mass = book<TH1F>("mass_"+name,"Mass "+label,50,minMass,maxMass);
  return hists;
}

template<typename T>
void BprimeHypHists::fill_BaseHists(const T & particle, BaseHists & hists, double weight){
  hists.pt->Fill(particle.pt(),weight);
  hists.eta->Fill(particle.eta(),weight);
  hists.phi->Fill(particle.phi(),weight);
  hists.mass->Fill(sqrt(particle.M2()),weight);
}

BprimeHypHists::BprimeHypHists(Context & ctx, const string & dirname, const string & hyp_name): Hists(ctx, dirname){
  wHad = book_BaseHists("wHad","W_{had} Hypothesis");
  wLep = book_BaseHists("wLep","W_{lep} Hypothesis"); 
  topLep = book_BaseHists("topLep","Top_{lep} Hypothesis");
  topHad = book_BaseHists("topHad","Top_{had} Hypothesis");
  mass = book_BaseHists("hyp","Hypothesis",50,3000); 
  mass_lep = book_BaseHists("Mass_lep","lep. Hypothesis",50,3000); 
  mass_had = book_BaseHists("Mass_had","had. Hypothesis",50,3000); 

  deltaR_w    = book<TH1F>("deltaR_w","#Delta R (W_{lep},W_{had})", 100, 0, 8);
  deltaPhi_w  = book<TH1F>("deltaPhi_w","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);
  deltaR_top  = book<TH1F>("deltaR_top","#Delta R (t_{lep},t_{had})", 100, 0, 4);
  deltaR_wtop = book<TH1F>("deltaR_wtop","#Delta R (t,W)", 100, 0, 5);

  pTratio_wtop   = book<TH1F>("pTratio_wtop","ratio pT W/Top", 100, 0, 10); 
  pTratio_toptop = book<TH1F>("pTratio_toptop","ratio pT Top_{lep}/Top_{had}", 100, 0, 10); 
  pTratio_ww     = book<TH1F>("pTratio_ww","ratio pT W_{lep}/W_{had}", 100, 0, 10); 

  recotype_h = book<TH1F>("recoType","0 hadronic Top, 1 leptonic Top",15, -0.5, 14.5);
  chiDis = book<TH1F>("chiVal","#Chi^{2} Value",100, 0,60);
  chiDis_lep = book<TH1F>("chiVal_lep","#Chi^{2} Value lep. Top",100, 0,60);
  chiDis_had = book<TH1F>("chiVal_had","#Chi^{2} Value had. Top",100, 0,60);

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

  chi_top_pT       = book<TH2F>("chi_top_pT","#Chi^{2} - Top p_{T}", 100, 0, 80, 100,0,1500);  
  chi_wlep_pT      = book<TH2F>("chi_wlep_pT","#Chi^{2} - W_{lep} p_{T}", 100, 0, 80, 100, 0,1500);  
  chi_whad_pT      = book<TH2F>("chi_whad_pT","#Chi^{2} - W_{had} p_{T}", 100, 0, 80, 100, 0,1500);  
  //chi_ST           = book<TH2F>("chi_ST","", 100, 0, 2, 100,0,800);  
  chi_deltaR_w_top = book<TH2F>("chi_deltaR_w_top","#Chi^{2} - #Delta R (t,W)", 100, 0, 80, 100,0,6);  
  


  recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
  gen = ctx.get_handle<BprimeGenContainer>("BprimeGen");
}

BprimeHypHists::~BprimeHypHists(){}


void BprimeHypHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  BprimeContainer hyp = event.get(recohyp);
  LorentzVector whad = hyp.get_wHad() ;
  LorentzVector wlep = hyp.get_wLep();
  LorentzVector topJets = hyp.get_topJets();
  LorentzVector thad = hyp.get_topHad();
  LorentzVector tlep = hyp.get_topLep();
  BprimeGenContainer GenInfo = event.get(gen);
  double chiVal = hyp.get_chiVal();
  int recotype = hyp.get_RecoTyp();

  //cout<<recotype<<endl;
  /*
  if(recotype==0){
    if((topJets_best+wlep_best).M()<100) return;
  }
  */
  //special treatment since topjets are not needed or filled for toptags
  if(recotype!=2)
    fill_BaseHists(topJets+whad+wlep, mass, weight);
  else if(recotype==2)
    fill_BaseHists(thad+wlep, mass, weight);
  
  fill_BaseHists(whad, wHad, weight);
  fill_BaseHists(wlep, wLep, weight);
  chiDis->Fill(chiVal,weight);
  deltaR_w->Fill(deltaR(whad,wlep),weight);
  deltaPhi_w->Fill(deltaPhi(whad,wlep),weight);
  recotype_h->Fill(recotype,weight);
  pTratio_ww->Fill(wlep.pt()/whad.pt(),weight);     
  chi_wlep_pT->Fill(chiVal,wlep.pt(),weight);
  chi_whad_pT->Fill(chiVal,whad.pt(),weight);
  //hadronic tops
  if(recotype==12 || recotype==0 || recotype==2){
    fill_BaseHists(thad,topHad,weight);
    fill_BaseHists(thad+wlep,mass_had,weight);
    deltaR_wtop->Fill(deltaR(thad,wlep),weight);
    chi_deltaR_w_top->Fill(chiVal,deltaR(thad,wlep),weight);
    pTratio_wtop->Fill(wlep.pt()/thad.pt(),weight);   
    chi_top_pT->Fill(chiVal,thad.pt(),weight);
    if(recotype==0) {
      deltaR_top->Fill(deltaR(thad,tlep),weight);
      pTratio_toptop->Fill(tlep.pt()/thad.pt(),weight);
    } 
    if(recotype==12) chiDis_had->Fill(chiVal,weight);
  }
  //leptonic tops
  if(recotype==11 || recotype==0){
    fill_BaseHists(tlep,topLep,weight);
    fill_BaseHists(tlep+whad,mass_lep,weight);
    deltaR_wtop->Fill(deltaR(tlep,whad),weight);
    pTratio_wtop->Fill(whad.pt()/tlep.pt(),weight);   
    chi_deltaR_w_top->Fill(chiVal,deltaR(tlep,whad),weight);
    chi_top_pT->Fill(chiVal,tlep.pt(),weight);
    if(recotype==11) chiDis_lep->Fill(chiVal,weight);
  }
  if(GenInfo.get_wHad().pt()>0){
    wHad_res_pt->Fill((GenInfo.get_wHad().pt()-whad.pt())/GenInfo.get_wHad().pt(),weight);   
    wHad_res_E->Fill((GenInfo.get_wHad().E()-whad.E())/GenInfo.get_wHad().E(),weight);       
    wHad_res_mass->Fill((GenInfo.get_wHad().M()-whad.M())/GenInfo.get_wHad().M(),weight);    
    wHad_res_phi->Fill(deltaPhi(GenInfo.get_wHad(),whad),weight);     
    wHad_res_eta->Fill(abs(GenInfo.get_wHad().eta()-whad.eta()),weight);     
    wHad_res_deltaR->Fill(deltaR(GenInfo.get_wHad(),whad),weight);
    wReco_dR_pTres_had->Fill(deltaR(GenInfo.get_wHad(),whad),(GenInfo.get_wHad().pt()-whad.pt())/GenInfo.get_wHad().pt(),weight);
    wReco_dR_pT_had->Fill(deltaR(GenInfo.get_wHad(),whad),GenInfo.get_wHad().pt(),weight);     
  }  
  if(GenInfo.get_wLep().pt()>0){
    wLep_res_pt->Fill((GenInfo.get_wLep().pt()-wlep.pt())/GenInfo.get_wLep().pt(),weight);     
    wLep_res_E->Fill((GenInfo.get_wLep().E()-wlep.E())/GenInfo.get_wLep().E(),weight);       
    wLep_res_mass->Fill((GenInfo.get_wLep().M()-wlep.M())/GenInfo.get_wLep().M(),weight);    
    wLep_res_phi->Fill(deltaPhi(GenInfo.get_wLep(),wlep),weight);     
    wLep_res_eta->Fill(abs(GenInfo.get_wLep().eta()-wlep.eta()),weight);     
    wLep_res_deltaR->Fill(deltaR(GenInfo.get_wLep(),wlep),weight);
    wReco_dR_pTres_lep->Fill(deltaR(GenInfo.get_wLep(),wlep),(GenInfo.get_wLep().pt()-wlep.pt())/GenInfo.get_wLep().pt(),weight);  
    wReco_dR_pT_lep->Fill(deltaR(GenInfo.get_wLep(),wlep),GenInfo.get_wLep().pt(),weight);       
  }
  if(GenInfo.get_topLep().pt()>0 && (recotype==11 || recotype==0)){
    topReco_dR_pT_lep->Fill(deltaR(GenInfo.get_topLep(),wlep+topJets),GenInfo.get_topLep().pt(),weight);   
    topReco_dR_pTres_lep->Fill(deltaR(GenInfo.get_topLep(),wlep+topJets),(GenInfo.get_topLep().pt()-(wlep+topJets).pt())/GenInfo.get_topLep().pt(),weight);
    topLep_res_pt->Fill((GenInfo.get_topLep().pt()-(wlep+topJets).pt())/GenInfo.get_topLep().pt(),weight);     
    topLep_res_E->Fill((GenInfo.get_topLep().E()-(wlep+topJets).E())/GenInfo.get_topLep().E(),weight);       
    topLep_res_mass->Fill((GenInfo.get_topLep().M()-(wlep+topJets).M())/GenInfo.get_topLep().M(),weight);    
    topLep_res_phi->Fill(deltaPhi(GenInfo.get_topLep(),wlep+topJets),weight);     
    topLep_res_eta->Fill(abs(GenInfo.get_topLep().eta()-(wlep+topJets).eta()),weight);     
    topLep_res_deltaR->Fill(deltaR(GenInfo.get_topLep(),wlep+topJets),weight);
  }
  if(GenInfo.get_topHad().pt()>0 && (recotype==12 || recotype==0 || recotype==2)){
    topReco_dR_pT_had->Fill(deltaR(GenInfo.get_topHad(),thad),GenInfo.get_topHad().pt(),weight);   
    topReco_dR_pTres_had->Fill(deltaR(GenInfo.get_topHad(),thad),(GenInfo.get_topHad().pt()-(thad).pt())/GenInfo.get_topHad().pt(),weight);
    topHad_res_pt->Fill((GenInfo.get_topHad().pt()-(thad).pt())/GenInfo.get_topHad().pt(),weight);   
    topHad_res_E->Fill((GenInfo.get_topHad().E()-(thad).E())/GenInfo.get_topHad().E(),weight);       
    topHad_res_mass->Fill((GenInfo.get_topHad().M()-(thad).M())/GenInfo.get_topHad().M(),weight);    
    topHad_res_phi->Fill(deltaPhi(GenInfo.get_topHad(),thad),weight);     
    topHad_res_eta->Fill(abs(GenInfo.get_topHad().eta()-(thad).eta()),weight);     
    topHad_res_deltaR->Fill(deltaR(GenInfo.get_topHad(),thad),weight);
  }
}


//double calc_Chi(LorentzVector top, LorentzVector w){
//}
