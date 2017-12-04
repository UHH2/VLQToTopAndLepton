#include "UHH2/VLQToTopAndLepton/include/BprimeHypHists.h"

using namespace std;
using namespace uhh2;

BprimeHypHists::BaseHists BprimeHypHists::book_BaseHists(const std::string & name, const std::string & label, double minMass, double maxMass, double minPt, double maxPt){
  BaseHists hists;
  hists.pt   = book<TH1F>("pt_"+name,"p_{T} "+ label+" [GeV]",100,minPt,maxPt);
  hists.eta  = book<TH1F>("eta_"+name,"#eta "+label,100,-4,4);
  hists.phi  = book<TH1F>("phi_"+name,"#phi "+label,100,-3.2,3.2);
  hists.mass = book<TH1F>("mass_"+name,"Mass "+label+" [GeV]",50,minMass,maxMass);
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
  //btagReco = book_BaseHists("btagReco","highest b-Tag",0,200);
  btagReco_csv = book<TH1F>("btagReco_cvs","highest b-Tag csv", 40, 0, 1);
  btagtoplepReco_csv = book<TH1F>("btagtoplepReco_cvs","highest b-Tag csv leptonic Top", 40, 0, 1);
  btagtophadReco_csv = book<TH1F>("btagtophadReco_cvs","highest b-Tag csv hadronic Top", 40, 0, 1);

  wHad = book_BaseHists("wHad","W_{had}",0,300);
  wLep = book_BaseHists("wLep","W_{lep}",0,300); 
  topLep = book_BaseHists("topLep","Top_{lep}",0,400);
  topHad = book_BaseHists("topHad","Top_{had}",0,400);
  mass = book_BaseHists("hyp","B",50,3000); 
  mass_lep = book_BaseHists("Mass_lep","lep. Hypothesis",50,3000); 
  mass_had = book_BaseHists("Mass_had","had. Hypothesis",50,3000); 

  forward = book_BaseHists("forward_jet","forward Jet",0,100,20,400);
  //balance = book_BaseHists("balance_jet","balance Jet",0,100,20,400); 
  combination = book_BaseHists("combination","forward + balance Jet",20,400,20,400); 
  no_forwardJet= book_BaseHists("no_forward_Jet","no forward Jet",50,3000); 

  deltaR_w    = book<TH1F>("deltaR_w","#Delta R (W_{lep},W_{had})", 100, 0, 8);
  deltaPhi_w  = book<TH1F>("deltaPhi_w","#Delta #phi (W_{lep},W_{had})", 100, 0, 4);
  deltaR_top  = book<TH1F>("deltaR_top","#Delta R (t_{lep},t_{had})", 100, 0, 4);
  deltaR_wtop = book<TH1F>("deltaR_wtop","#Delta R (t,W)", 100, 0, 5);

  pTratio_wtop   = book<TH1F>("pTratio_wtop","ratio pT W/Top", 100, 0, 10); 
  pTratio_toptop = book<TH1F>("pTratio_toptop","ratio pT Top_{lep}/Top_{had}", 100, 0, 10); 
  pTratio_ww     = book<TH1F>("pTratio_ww","ratio pT W_{lep}/W_{had}", 100, 0, 10); 

  recotype_h = book<TH1F>("recoType","reconstruction typ",15, -0.5, 14.5);
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

  Bprime_res_pt   = book<TH1F>("Bprime_res_pt","B Resolution p_{T}", 100, -3, 2);
  Bprime_res_E    = book<TH1F>("Bprime_res_E","B Resolution E", 100, -3, 2); 
  Bprime_res_mass = book<TH1F>("Bprime_res_mass","B Resolution Mass", 100, -3, 2); 
  Bprime_res_phi  = book<TH1F>("Bprime_res_phi","B Resolution #phi", 100, 0, 3.14); 
  Bprime_res_eta  = book<TH1F>(" Bprime_res_eta","B Resolution #eta", 100, 0, 4); 
  Bprime_res_deltaR  = book<TH1F>(" Bprime_res_deltaR","B Resolution #Delta R", 100, 0, 4); 

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
  
  //forward Jet TH1F/TH2F
  deltaR_forward_B = book<TH1F>("deltaR_forward_B","#Delta R (Jet_{forward},whad)", 100, 0,8); 
  forward_pt_eta = book<TH2F>("forward_pt_eta","Jet_{forward} pT #eta", 100, 0, 250, 100,-8,8);  
  
  //matching studies
  matched_top_lep          = book_BaseHists("matched_top_lep","Top_{lep,match}");
  matched_top_had          = book_BaseHists("matched_top_had","Top_{had,match}");
  matched_W_lep            = book_BaseHists("matched_W_lep","W_{lep,match}");
  matched_W_had            = book_BaseHists("matched_W_had","W_{had,match}");
  matched_tops             = book_BaseHists("matched_tops","matched Tops");
  semi_matched_tops        = book_BaseHists("semi_matched_tops","semi matched Tops");
  matched_top_matched_W    = book_BaseHists("matched_top_matched_W","Top_{match} with W_{match}");
  matched_W_matched_top    = book_BaseHists("matched_W_matched_top","W_{match} with Top{match}");
  matched_top_unmatched_W  = book_BaseHists("matched_top_unmatched_W","Top_{match} with W_{unmatch}");
  matched_W_unmatched_top  = book_BaseHists("matched_W_unmatched_top","W_{match} with Top_{unmatch}");

  reco_pt_tW    = book<TH2F>("reco_pt_tW","Reco. p_{T,t} p_{T,W}", 100, 0, 1000,100,0,1000);
  gen_pt_tW    = book<TH2F>("gen_pt_tW","Gen. p_{T,t} p_{T,W}", 100, 0, 1000,100,0,1000);

  deltaPhi_tlep_whad = book<TH1F>("deltaPhi_tlep_whad","#Delta #phi (t_{lep},W_{had})", 100, 0, 8);
  deltaEta_tlep_whad = book<TH1F>("deltaEta_tlep_whad","#delta #eta (t_{lep},W_{had})", 100, 0, 8);
  deltaPhi_thad_wlep = book<TH1F>("deltaPhi_thad_wlep","#Delta #phi (t_{had},W_{lep})", 100, 0, 8);
  deltaEta_thad_wlep = book<TH1F>("deltaEta_thad_wlep","#Delta #eta (t_{had},W_{lep})", 100, 0, 8);

  deltaPhi_deltaEta_tlep_whad = book<TH2F>("deltaPhi_deltaEta_tlep_whad","", 100, 0, 8, 100,0,8); 
  deltaPhi_deltaEta_thad_wlep = book<TH2F>("deltaPhi_deltaEta_thad_wlep","", 100, 0, 8, 100,0,8); 
  
  
  dR_forwardJet_bprime = book<TH1F>("dR_forwardJet_bprime","#Delta R (forwardJet,B)", 100, 0, 8);
  deltaPhi_forwardJet_bprime = book<TH1F>("deltaPhi_forwardJet_bprime","#Delta #phi (forward Jet,B)", 100, 0, 3.4);
  deltaEta_forwardJet_bprime = book<TH1F>("deltaEta_forwardJet_bprime","#Delta #eta (forward Jet,B)", 100, 0, 8);
  dRmin_forwardJet_top_w = book<TH1F>("dRmin_forwardJet_top_w","#Delta Rmin(forward Jet,(top,W))", 100, 0, 8);
  dR_forwardJet_top = book<TH1F>("dR_forwardJet_top","#Delta R (forward Jet,top)", 100, 0, 8);
  dR_forwardJet_w = book<TH1F>("dR_forwardJet_w","#Delta R (forward Jet,W)", 100, 0, 8);
    
  deltaEta_bprime_forwardJet_eta = book<TH2F>("deltaEta_bprime_forwardJet_eta","#Delta #eta (forward Jet,B) forward Jet #eta", 100, 0, 8, 100,-5,5);

  forwardJet_bprime_deltaPhi_deltaEta = book<TH2F>("forwardJet_bprime_deltaPhi_deltaEta","#Delta #phi #Delta #eta (forwardJet,B)", 100, 0, 3.4, 100, 0, 8);

  recohyp = ctx.get_handle<BprimeContainer>(hyp_name);
  gen = ctx.get_handle<BprimeGenContainer>("BprimeGen");
  /*
  weight_toptag = ctx.get_handle<double>("weight_toptag");
  weight_toptag_up = ctx.get_handle<double>("weight_toptag_up");
  weight_toptag_down = ctx.get_handle<double>("weight_toptag_up");
  */
}

BprimeHypHists::~BprimeHypHists(){}

void BprimeHypHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  BprimeContainer hyp = event.get(recohyp);
  LorentzVector whad = hyp.get_wHad();
  LorentzVector wlep = hyp.get_wLep();
  LorentzVector topJets = hyp.get_topJets();
  LorentzVector thad = hyp.get_topHad();
  LorentzVector tlep = hyp.get_topLep();
  LorentzVector bprime;
  std::vector<LorentzVector> jets_top = hyp.get_topLorentz();
  std::vector<LorentzVector> jets_whad = hyp.get_wHadLorentz();
  double chiVal = hyp.get_chiVal();
  int recotype = hyp.get_RecoTyp();
  //vector<Jet>* jets = event.jets;

  double matching_distance = 0.4;

  double btagdis = hyp.get_btag_discriminator();
  double btagwhad = hyp.get_btag_whad();
  double toptag_weight = -1;

  if(btagdis != -1)
    btagReco_csv->Fill(btagdis,weight);
  /*
  if(recotype ==2 && !event.isRealData){
    if(thad.pt()>=400 && thad.pt()<=550){
      //event.weight *= 0.91;
      weight *= 0.91;
      toptag_weight = 0.91;
    }
    else if(thad.pt()>550){
      //event.weight *= 0.97;
      weight *= 0.97;
      toptag_weight = 0.97;
    }
    }*/

  /*
  event.set(weight_toptag,toptag_weight);
  event.set(weight_toptag_up,toptag_weight*toptag_weight);
  event.set(weight_toptag_down,event.weight);
  */
  /*
  if(whad.M2()>30)
    whad.SetE(sqrt(80.4*80.4+whad.Px()*whad.Px()+whad.Py()*whad.Py()+whad.Pz()*whad.Pz()));
  if(thad.M2()>30)
    thad.SetE(sqrt(175*175+thad.Px()*thad.Px()+thad.Py()*thad.Py()+thad.Pz()*thad.Pz()));
  if(tlep.M2()>30)
    tlep.SetE(sqrt(175*175+tlep.Px()*tlep.Px()+tlep.Py()*tlep.Py()+tlep.Pz()*tlep.Pz()));
  */
 

  //cout<<recotype<<endl;
  /*
  if(recotype==0){
    if((topJets_best+wlep_best).M()<100) return;
  }
  */
  //special treatment since topjets are not needed or filled for toptags

  if(recotype==11||recotype==6){
    bprime = tlep+whad;
    btagtoplepReco_csv->Fill(btagwhad,weight);
    reco_pt_tW->Fill(tlep.pt(),whad.pt(),weight);
    //cout<< "tlep,Whad delta phi "<<deltaPhi(tlep,whad)<<" delta eta "<<abs(tlep.eta()-whad.eta())<<endl; 
    deltaPhi_tlep_whad->Fill(deltaPhi(tlep,whad),weight);
    deltaEta_tlep_whad->Fill(abs(tlep.eta()-whad.eta()),weight);
    deltaPhi_deltaEta_tlep_whad->Fill(deltaPhi(tlep,whad),abs(tlep.eta()-whad.eta()),weight);
  }
  else if(recotype==12 || recotype==2){
    bprime = thad+wlep;
    btagtophadReco_csv->Fill(btagwhad,weight);
    if(btagdis==-1){
      btagReco_csv->Fill(btagwhad,weight);
    }
    reco_pt_tW->Fill(thad.pt(),wlep.pt(),weight);
    deltaPhi_thad_wlep->Fill(deltaPhi(thad,wlep),weight);
    deltaEta_thad_wlep->Fill(abs(thad.eta()-wlep.eta()),weight);
    deltaPhi_deltaEta_thad_wlep->Fill(deltaPhi(thad,wlep),abs(thad.eta()-wlep.eta()),weight);
    //cout<< "thad,Wlep delta phi "<<deltaPhi(thad,wlep)<<" delta eta "<<abs(thad.eta()-wlep.eta())<<endl; 
  }
  
  if(recotype!=2)
    fill_BaseHists(topJets+whad+wlep, mass, weight);
  else if(recotype==2)
    fill_BaseHists(thad+wlep, mass, weight);

  if((btagdis==-1 && recotype==12 && btagwhad==-1) || (btagdis==-1 && recotype!=12))
    btagReco_csv->Fill(btagwhad,weight);
  

  if(jets_top.size() ==1 && jets_whad.size() ==2 && 1==2){
    cout<<"======================================================================================="<<endl;
    cout<<"Reco Type "<<recotype<<" chi2 "<<chiVal<<" B Mass "<<bprime.M()<<" Top lep Mass "<<tlep.M()<<" W had Mass "<<whad.M()<<endl;
    cout<<"W_lep: Mass "<<wlep.M()<<" ("<<wlep.px()<<","<<wlep.py()<<","<<wlep.pz()<<","<<wlep.energy()<<")"<<endl; 
    cout<<"Top_jets ("<<jets_top[0].px()<<","<<jets_top[0].py()<<","<<jets_top[0].pz()<<","<<jets_top[0].energy()<<")"<<endl; 
    cout<<"W_had :"<<endl;
    cout<<"jet1 ("<<jets_whad[0].px()<<","<<jets_whad[0].py()<<","<<jets_whad[0].pz()<<","<<jets_whad[0].energy()<<")"<<endl; 
    cout<<"jet2 ("<<jets_whad[1].px()<<","<<jets_whad[1].py()<<","<<jets_whad[1].pz()<<","<<jets_whad[1].energy()<<")"<<endl; 
  }

  //if( recotype == 12 || recotype == 11)
  //cout<<"recotype "<<recotype<< " "<<bprime.M()<<" "<<hyp.get_Mass()<<endl;



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
  if(recotype==11 || recotype==0 || recotype==6){
    fill_BaseHists(tlep,topLep,weight);
    fill_BaseHists(tlep+whad,mass_lep,weight);
    deltaR_wtop->Fill(deltaR(tlep,whad),weight);
    pTratio_wtop->Fill(whad.pt()/tlep.pt(),weight);   
    chi_deltaR_w_top->Fill(chiVal,deltaR(tlep,whad),weight);
    chi_top_pT->Fill(chiVal,tlep.pt(),weight);
    if(recotype==11) chiDis_lep->Fill(chiVal,weight);
  }

  //most forward jet
  LorentzVector forwardJet (0,0,0,0);
  LorentzVector balanceJet (0,0,0,0);
  for(auto jet : *event.jets)
    if(abs(forwardJet.eta())<abs(jet.eta()))
      forwardJet = jet.v4();
  /*
  //forward Jet
  LorentzVector forwardJet (0,0,0,0);
  LorentzVector balanceJet (0,0,0,0);
  if(recotype==11 || recotype==12||recotype==0){
    string unusedJets = hyp.get_unusedJets();
    //cout<<unusedJets.size()<<" "<<jets->size()<<endl;
    for(unsigned int i =0; i<jets->size();i++){ 
      if(!unusedJets[i]){
	if(fabs(forwardJet.eta())<fabs(jets->at(i).eta())){
	  balanceJet = forwardJet;
	  forwardJet = jets->at(i).v4();
	}
      }
    }
    //cout<<"forwardJet eta "<< forwardJet.eta()<<" balanceJet eta "<<balanceJet.eta()<<endl;
  }
  else if(recotype==2 || recotype==6){
    for(unsigned int i =0; i<jets->size();i++){ 
      if(deltaR(jets->at(i).v4(),thad)>2)
	if(fabs(forwardJet.eta())<fabs(jets->at(i).eta())){
	  balanceJet = forwardJet;
	  forwardJet = jets->at(i).v4();
	}
    }
  }

  */
  
  if(forwardJet.pt()>0){
    if(recotype == 11 || recotype==0 || recotype==6){
      dR_forwardJet_top->Fill(deltaR(forwardJet,tlep),weight);
      dR_forwardJet_w->Fill(deltaR(forwardJet,whad),weight);
      dRmin_forwardJet_top_w->Fill(deltaR(forwardJet,tlep) < deltaR(forwardJet,whad) ? deltaR(forwardJet,tlep) : deltaR(forwardJet,whad), weight);
    }
    if(recotype==12 || recotype==0 || recotype==2){
      dR_forwardJet_top->Fill(deltaR(forwardJet,thad),weight);
      dR_forwardJet_w->Fill(deltaR(forwardJet,wlep),weight);
      dRmin_forwardJet_top_w->Fill(deltaR(forwardJet,thad) < deltaR(forwardJet,wlep) ? deltaR(forwardJet,thad) : deltaR(forwardJet,wlep), weight);				       
    }
    dR_forwardJet_bprime->Fill(deltaR(bprime,forwardJet),weight);
    deltaPhi_forwardJet_bprime->Fill(deltaPhi(bprime,forwardJet),weight);
    deltaEta_forwardJet_bprime->Fill(abs(bprime.eta()-forwardJet.eta()),weight);
    forwardJet_bprime_deltaPhi_deltaEta->Fill(deltaPhi(bprime,forwardJet),abs(bprime.eta()-forwardJet.eta()),weight);
    fill_BaseHists(forwardJet,forward,weight);
    //fill_BaseHists(forwardJet+balanceJet+whad+wlep+topJets,combination,weight);
    deltaR_forward_B->Fill(deltaR(forwardJet,whad),weight);
    forward_pt_eta->Fill(forwardJet.pt(),forwardJet.eta(),weight);
    deltaEta_bprime_forwardJet_eta->Fill(deltaPhi(bprime,forwardJet),bprime.eta(),weight);
  }
  //else
  //  fill_BaseHists(whad+wlep+topJets,no_forwardJet,weight);
  

  /*
  if(balanceJet.pt()>0)fill_BaseHists(balanceJet,balance,weight);
  */

  if(event.isRealData) return;
  BprimeGenContainer GenInfo = event.get(gen);
  //cout<<" Reco Type "<<recotype << " pT B "<<bprime.pt()<< " t lep "<<tlep.pt()<<" w had "<< whad.pt()<< " t had "<<thad.pt()<<" w lep "<< wlep.pt()<<" Gen t lep "<<GenInfo.get_topLep().pt()<<" w had "<< GenInfo.get_wHad().pt()<<" t had "<<GenInfo.get_topHad().pt()<<" w lep "<< GenInfo.get_wLep().pt() <<endl;
  
  if(GenInfo.get_wHad().pt()>0){
    wHad_res_pt->Fill((GenInfo.get_wHad().pt()-whad.pt())/GenInfo.get_wHad().pt(),weight);   
    wHad_res_E->Fill((GenInfo.get_wHad().E()-whad.E())/GenInfo.get_wHad().E(),weight);       
    wHad_res_mass->Fill((GenInfo.get_wHad().M()-whad.M())/GenInfo.get_wHad().M(),weight);    
    wHad_res_phi->Fill(deltaPhi(GenInfo.get_wHad(),whad),weight);     
    wHad_res_eta->Fill(abs(GenInfo.get_wHad().eta()-whad.eta()),weight);     
    wHad_res_deltaR->Fill(deltaR(GenInfo.get_wHad(),whad),weight);
    wReco_dR_pTres_had->Fill(deltaR(GenInfo.get_wHad(),whad),(GenInfo.get_wHad().pt()-whad.pt())/GenInfo.get_wHad().pt(),weight);
    wReco_dR_pT_had->Fill(deltaR(GenInfo.get_wHad(),whad),GenInfo.get_wHad().pt(),weight); 

    if(deltaR(GenInfo.get_wHad(), whad)<matching_distance && whad.pt()>0 && GenInfo.get_wHad().pt()>0){
      fill_BaseHists(whad,matched_W_had, weight);
      if(tlep.pt()>0 && deltaR(tlep,GenInfo.get_topLep())>= 0.4 && GenInfo.get_topLep().pt()>0)
	fill_BaseHists(whad,matched_W_unmatched_top,weight);
    }  
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

    if(deltaR(GenInfo.get_wLep(), wlep)<matching_distance && wlep.pt()>0 && GenInfo.get_wLep().pt()>0){
      fill_BaseHists(wlep,matched_W_lep, weight);
      if(thad.pt()>0 && deltaR(thad,GenInfo.get_topHad())>= 0.4 && GenInfo.get_topHad().pt()>0)
	fill_BaseHists(wlep,matched_W_unmatched_top,weight);
    }

  }
  if(GenInfo.get_topLep().pt()>0 && (recotype==11 || recotype==0 || recotype==6)){
    gen_pt_tW->Fill(GenInfo.get_topLep().pt(),GenInfo.get_wHad().pt(),weight);
    topReco_dR_pT_lep->Fill(deltaR(GenInfo.get_topLep(),wlep+topJets),GenInfo.get_topLep().pt(),weight);   
    topReco_dR_pTres_lep->Fill(deltaR(GenInfo.get_topLep(),wlep+topJets),(GenInfo.get_topLep().pt()-(wlep+topJets).pt())/GenInfo.get_topLep().pt(),weight);
    topLep_res_pt->Fill((GenInfo.get_topLep().pt()-(wlep+topJets).pt())/GenInfo.get_topLep().pt(),weight);     
    topLep_res_E->Fill((GenInfo.get_topLep().E()-(wlep+topJets).E())/GenInfo.get_topLep().E(),weight);       
    topLep_res_mass->Fill((GenInfo.get_topLep().M()-(wlep+topJets).M())/GenInfo.get_topLep().M(),weight);    
    topLep_res_phi->Fill(deltaPhi(GenInfo.get_topLep(),wlep+topJets),weight);     
    topLep_res_eta->Fill(abs(GenInfo.get_topLep().eta()-(wlep+topJets).eta()),weight);     
    topLep_res_deltaR->Fill(deltaR(GenInfo.get_topLep(),wlep+topJets),weight);

    if(deltaR(GenInfo.get_topLep(),tlep)<matching_distance && tlep.pt()>0){
      fill_BaseHists(tlep,matched_top_lep, weight);
      if(GenInfo.get_topHad().pt()>0 && deltaR(GenInfo.get_topHad(),thad)<matching_distance && thad.pt()>0){
	fill_BaseHists(tlep,matched_tops, weight);
	fill_BaseHists(thad,matched_tops, weight);
      }
      else if(GenInfo.get_topHad().pt()>0 && thad.pt()>0 )
	fill_BaseHists(tlep,semi_matched_tops,weight);
      if(GenInfo.get_wHad().pt()>0 && whad.pt()>0){
	if(deltaR(GenInfo.get_wHad(),whad)<matching_distance){
	  fill_BaseHists(thad,matched_top_matched_W, weight);
	  fill_BaseHists(whad,matched_W_matched_top, weight);
	}
	else{
	  fill_BaseHists(thad,matched_top_unmatched_W,weight); 
	}
      }
    }
  }
  if(GenInfo.get_topHad().pt()>0 && (recotype==12 || recotype==0 || recotype==2)){
    gen_pt_tW->Fill(GenInfo.get_topHad().pt(),GenInfo.get_wLep().pt(),weight);
    topReco_dR_pT_had->Fill(deltaR(GenInfo.get_topHad(),thad),GenInfo.get_topHad().pt(),weight);   
    topReco_dR_pTres_had->Fill(deltaR(GenInfo.get_topHad(),thad),(GenInfo.get_topHad().pt()-(thad).pt())/GenInfo.get_topHad().pt(),weight);
    topHad_res_pt->Fill((GenInfo.get_topHad().pt()-(thad).pt())/GenInfo.get_topHad().pt(),weight);   
    topHad_res_E->Fill((GenInfo.get_topHad().E()-(thad).E())/GenInfo.get_topHad().E(),weight);       
    topHad_res_mass->Fill((GenInfo.get_topHad().M()-(thad).M())/GenInfo.get_topHad().M(),weight);    
    topHad_res_phi->Fill(deltaPhi(GenInfo.get_topHad(),thad),weight);     
    topHad_res_eta->Fill(abs(GenInfo.get_topHad().eta()-(thad).eta()),weight);     
    topHad_res_deltaR->Fill(deltaR(GenInfo.get_topHad(),thad),weight);

    if(deltaR(GenInfo.get_topHad(), thad)<matching_distance){
      fill_BaseHists(thad,matched_top_had, weight);
      if(GenInfo.get_topLep().pt()>0 && deltaR(GenInfo.get_topLep(),thad)>=matching_distance && tlep.pt()>0)
	 fill_BaseHists(thad,semi_matched_tops,weight);
      if(GenInfo.get_wLep().pt()>0 && wlep.pt()>0){
	if(deltaR(GenInfo.get_wLep(),wlep)<matching_distance){
	  fill_BaseHists(thad,matched_top_matched_W, weight);
	  fill_BaseHists(wlep,matched_W_matched_top, weight);
	}
	else{
	  fill_BaseHists(thad,matched_top_unmatched_W,weight); 
	}
      }
    }
  }
  if(GenInfo.get_bprime().pt()>0){
    Bprime_res_pt->Fill((GenInfo.get_bprime().pt()-bprime.pt())/GenInfo.get_bprime().pt(),weight);   
    Bprime_res_E->Fill((GenInfo.get_bprime().E()-bprime.E())/GenInfo.get_bprime().E(),weight);   
    Bprime_res_mass->Fill((GenInfo.get_bprime().M()-bprime.M())/GenInfo.get_bprime().M(),weight); 
    Bprime_res_phi->Fill(deltaPhi(GenInfo.get_bprime(),bprime),weight);  
    Bprime_res_eta->Fill(abs(GenInfo.get_bprime().eta()-bprime.eta()),weight);  
    Bprime_res_deltaR->Fill(deltaR(GenInfo.get_bprime(),bprime),weight);  
   }

}

//double calc_Chi(LorentzVector top, LorentzVector w){
//}
