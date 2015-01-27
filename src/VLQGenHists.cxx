//framework includes
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"
#include "UHH2/core/include/Event.h"

//general includes
#include <iostream>

using namespace std;
using namespace uhh2;


VLQGenHists::VLQGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  
  //histograms for vector like quarks
  VLQ_Number =  book<TH1F>("VLQ_number", "Number of VLQ", 5, -0.5, 4.5);
  VLQ_pdgId =  book<TH1F>("VLQ_pdgId", "pdgId", 61, -30.5, 30.5); 

  VLQ_eta_lead = book<TH1F>("VLQ_eta_lead", "#eta_{VLQ}(lead)", 100, -4, 4);
  VLQ_eta_subl = book<TH1F>("VLQ_eta_subl", "#eta_{VLQ}(sublead)", 100, -4, 4);
  
  VLQ_phi_lead = book<TH1F>("VLQ_phi_lead", "#phi_{VLQ}(lead)", 64, -3.2, 3.2);
  VLQ_phi_subl = book<TH1F>("VLQ_phi_subl", "#phi_{VLQ}(sublead)", 64, -3.2, 3.2);
  
  VLQ_pt_lead  =  book<TH1F>("VLQ_pt_lead" , "p_{T}^{VLQ}(lead) [GeV/c]", 200, 0, 2000);
  VLQ_pt_subl  =  book<TH1F>("VLQ_pt_subl" , "p_{T}^{VLQ}(sublead) [GeV/c]", 200, 0, 2000);
  
  VLQ_decay    = book<TH1F>("VLQ_decay", "VLQ decay modes", 61, -30.5, 30.5);


  //number of bosons 
  NHiggs = book<TH1F>("NHiggs" , "Number of Higgs", 5, -0.5, 4.5);
  NW     = book<TH1F>("NW"     , "Number of W", 5, -0.5, 4.5);
  NZ     = book<TH1F>("NZ"     , "Number of Z", 5, -0.5, 4.5);

  //number of quarks
  Nbottom = book<TH1F>("Nbottom", "Number of b", 9, -.5, 8.5);
  Ntop    = book<TH1F>("Ntop"   , "Number of top", 6, -0.5, 5.5);

  //number of leptons
  Nlept      = book<TH1F>("Nlept"     , "Number of leptons", 7, -.5, 6.5);
  Nmu        = book<TH1F>("Nmu"       , "Number of muons", 5, -.5, 4.5);
  Nelectrons = book<TH1F>("Nelectrons", "Number of electrons",5, -.5, 4.5);
  
  //higgs 
  higgs_decay    = book<TH1F>("higgs_decay"     , "Higgs decay modes", 61, -30.5, 30.5);

  higgs_pt       = book<TH1F>("higgs_pt" , "p_{T}^{higgs} [GeV/c]", 200, 0, 2000);
  higgs_eta      = book<TH1F>("higgs_eta", "#eta_{higgs}", 40, -2.5, 2.5);
  higgs_phi      = book<TH1F>("higgs_phi", "#phi_{higgs}", 64, -3.2, 3.2);
    
  higgs_pt_lead  = book<TH1F>("higgs_pt_lead" , "p_{T}^{higgs}(lead) [GeV/c]", 200, 0, 2000);
  higgs_pt_subl  = book<TH1F>("higgs_pt_subl" , "p_{T}^{higgs}(sublead) [GeV/c]", 200, 0, 2000);
  higgs_eta_lead = book<TH1F>("higgs_eta_lead", "#eta_{higgs}(lead)", 40, -2.5, 2.5);
  higgs_eta_subl = book<TH1F>("higgs_eta_subl", "#eta_{higgs}(sublead)", 40, -2.5, 2.5);
  higgs_phi_lead = book<TH1F>("higgs_phi_lead", "#phi_{higgs}(lead)", 64, -3.2, 3.2);
  higgs_phi_subl = book<TH1F>("higgs_phi_subl", "#phi_{higgs}(sublead)", 64, -3.2, 3.2);

  DeltaR_bb      = book<TH1F>("DeltaR_bb"   , "#Delta R_{bb}", 50, 0, 5); //b variable

  //W
  W_decay    = book<TH1F>("W_decay"   , "W decay modes", 61, -30.5, 30.5);
  
  W_pt       = book<TH1F>("W_pt" , "p_{T}^{W} [GeV/c]", 200, 0, 2000);
  W_eta      = book<TH1F>("W_eta", "#eta_{W}", 40, -2.5, 2.5);     
  W_phi      = book<TH1F>("W_phi", "#phi_{W}", 64, -3.2, 3.2);
   
  W_pt_lead  = book<TH1F>("W_pt_lead" , "p_{T}^{W}(lead) [GeV/c]", 200, 0, 2000);
  W_pt_subl  = book<TH1F>("W_pt_subl" , "p_{T}^{W}(sublead) [GeV/c]", 200, 0, 2000);
  W_eta_lead = book<TH1F>("W_eta_lead", "#eta_{W}(lead)", 40, -2.5, 2.5);
  W_eta_subl = book<TH1F>("W_eta_subl", "#eta_{W}(sublead)", 40, -2.5, 2.5);
  W_phi_lead = book<TH1F>("W_phi_lead", "#phi_{W}(lead)", 64, -3.2, 3.2);
  W_phi_subl = book<TH1F>("W_phi_subl", "#phi_{W}(sublead)", 64, -3.2, 3.2);

  //Z 
  Z_decay    = book<TH1F>("Z_decay"   , "Z decay modes", 61, -30.5, 30.5);
  
  Z_pt       = book<TH1F>("Z_pt" , "p_{T}^{Z} [GeV/c]", 200, 0, 2000);
  Z_eta      = book<TH1F>("Z_eta", "#eta_{Z}", 40, -2.5, 2.5);
  Z_phi      = book<TH1F>("Z_phi", "#phi_{Z}", 64, -3.2, 3.2);
  
  Z_pt_lead  = book<TH1F>("Z_pt_lead" , "p_{T}^{Z}(lead) [GeV/c]", 200, 0, 2000);
  Z_pt_subl  = book<TH1F>("Z_pt_subl" , "p_{T}^{Z}(sublead) [GeV/c]", 200, 0, 2000);
  Z_eta_lead = book<TH1F>("Z_eta_lead", "#eta_{Z}(lead)", 40, -2.5, 2.5);
  Z_eta_subl = book<TH1F>("Z_eta_subl", "#eta_{Z}(sublead)", 40, -2.5, 2.5);
  Z_phi_lead = book<TH1F>("Z_phi_lead", "#phi_{Z}(lead)", 64, -3.2, 3.2);
  Z_phi_subl = book<TH1F>("Z_phi_subl", "#phi_{Z}(sublead)", 64, -3.2, 3.2);


  //top 
  top_decay    = book<TH1F>("top_decay"   , "Top decay modes", 61, -30.5, 30.5);
  
  top_pt       = book<TH1F>("top_pt" , "p_{T}^{top} [GeV/c]", 200, 0, 2000);
  top_eta      = book<TH1F>("top_eta", "#eta_{top}(lead)", 40, -2.5, 2.5);
  top_phi      = book<TH1F>("top_phi", "#phi_{top}(lead)", 64, -3.2, 3.2);
  
  top_pt_lead  = book<TH1F>("top_pt_lead" , "p_{T}^{top}(lead) [GeV/c]", 200, 0, 2000);
  top_pt_subl  = book<TH1F>("top_pt_subl" , "p_{T}^{top}(sublead) [GeV/c]", 200, 0, 2000);
  top_eta_lead = book<TH1F>("top_eta_lead", "#eta_{top}(lead)", 40, -2.5, 2.5);
  top_eta_subl = book<TH1F>("top_eta_subl", "#eta_{top}(sublead)", 40, -2.5, 2.5);
  top_phi_lead = book<TH1F>("top_phi_lead", "#phi_{top}(lead)", 64, -3.2, 3.2);
  top_phi_subl = book<TH1F>("top_phi_subl", "#phi_{top}(sublead)", 64, -3.2, 3.2);
  


  //no mothers
  particles_noMother     = book<TH1F>("particle_noMother"   , "Particles with no mother", 61, -30.5, 30.5);
  particles_noMother_pT  = book<TH1F>("particles_noMother_pT" , "p_{T} [GeV/c]", 200, 0, 1000);
  particles_noMother_eta = book<TH1F>("particles_noMother_eta", "#eta ", 40, -2.5, 2.5);
  particles_noMother_phi = book<TH1F>("particles_noMother_phi", "#phi", 64, -3.2, 3.2);

  //no daughters or mothers
  particles_noMotherNoDaughter     = book<TH1F>("particle_noMotherNoDaughter"   , "Particles with no mother or daughter", 61, -30.5, 30.5);
  particles_noMotherNoDaughter_pT  = book<TH1F>("particles_noMotherNoDaughter_pT" , "p_{T} [GeV/c]", 200, 0, 1000);
  particles_noMotherNoDaughter_eta = book<TH1F>("particles_noMotherNoDaughter_eta", "#eta ", 40, -2.5, 2.5);
  particles_noMotherNoDaughter_phi = book<TH1F>("particles_noMotherNoDaughter_phi", "#phi", 64, -3.2, 3.2);

}


void VLQGenHists::fill(const Event & event){
  // fill the histograms. Don't forget to always use the weight when filling:
  //     double weight = event.weight;
  double weight = event.weight;
  
  // for (unsigned int i=0; i<event.get_current_triggernames().size();i++)
  //  cout<< event.get_current_triggernames()[i]<<endl;

     
  const vector<GenParticle> * genparticles = event.genparticles;
  
  vector<GenParticle> vlq;
  
  vector<GenParticle> higgs;
  vector<GenParticle> wboson;
  vector<GenParticle> zboson;
  
  vector<GenParticle> top;
  vector<GenParticle> bquarks;
  vector<GenParticle> electrons;
  vector<GenParticle> muons;
  
  
  for(auto igenp : *genparticles){ 
      

    //cout<<"pdgId GenParticle: "<<igenp.pdgId()<<endl;


    //put all the particles in vectors
    if (abs(igenp.pdgId()) == 6000008 || abs(igenp.pdgId()) == 6000007)
      vlq.push_back(igenp);
    else if(abs(igenp.pdgId()) == 5 )
      bquarks.push_back(igenp);
    else if(abs(igenp.pdgId()) == 6 )
      top.push_back(igenp);
    else if (abs(igenp.pdgId()) == 11 )
      electrons.push_back(igenp);
    else if (abs(igenp.pdgId()) == 13 )
      muons.push_back(igenp);
    else if(abs(igenp.pdgId()) == 23 )
      zboson.push_back(igenp);
    else if(abs(igenp.pdgId()) == 24 )
      wboson.push_back(igenp);
    else if(abs(igenp.pdgId()) == 25 )
      higgs.push_back(igenp);
  
    
    const GenParticle* mother1 = igenp.mother(genparticles,1);
    const GenParticle* mother2 = igenp.mother(genparticles,2);

    int mother1_pdgId = 0; 
    if(mother1)mother1_pdgId = mother1->pdgId();
    int mother2_pdgId = 0;
    if(mother2) mother2_pdgId = mother2->pdgId();

    
    const GenParticle* daughter1 = igenp.daughter(genparticles,1);
    const GenParticle* daughter2 = igenp.daughter(genparticles,2);

    int daughter1_pdgId = 0; 
    if(daughter1)daughter1_pdgId = daughter1->pdgId();
    int daughter2_pdgId = 0;
    if(daughter2) daughter2_pdgId = daughter2->pdgId();


    // if(abs(igenp.pdgId())==25 ||abs(igenp.pdgId())==23)
    // cout<<"pdgId GenParticle: "<<igenp.pdgId() <<" mom1: "<<  mother1_pdgId<<" mom2: "<< mother2_pdgId <<" daughter1: "<<  daughter1_pdgId<<" daughter2: "<< daughter2_pdgId<<endl;
  
    //Fill decay histograms
    if(abs(mother1_pdgId)== 6000007 || abs(mother1_pdgId)== 6000008 || abs(mother2_pdgId)== 6000008 || abs(mother2_pdgId)== 6000007 ){
      VLQ_decay->Fill(igenp.pdgId());
      //cout<<igenp.pdgId()<<endl;
    }
    if(abs(mother1_pdgId)== 6 || abs(mother2_pdgId)== 6)
      top_decay->Fill(igenp.pdgId()); 
    if(abs(mother1_pdgId)== 23 || abs(mother2_pdgId)== 23)
      Z_decay->Fill(igenp.pdgId()); 
    if(abs(mother1_pdgId)== 24 || abs(mother2_pdgId)== 24)
      W_decay ->Fill(igenp.pdgId()); 
    if(abs(mother1_pdgId)== 25 || abs(mother2_pdgId)== 25)
      higgs_decay->Fill(igenp.pdgId());


    if(abs(mother1_pdgId)== 0 && abs(mother2_pdgId)== 0){
      particles_noMother->Fill(igenp.pdgId());           
      particles_noMother_pT->Fill(igenp.pt());  
      particles_noMother_eta->Fill(igenp.eta());  
      particles_noMother_phi->Fill(igenp.phi());   
    }
    if(abs(mother1_pdgId)== 0 && abs(mother2_pdgId)== 0 &&  abs(daughter2_pdgId)== 0 &&  abs(daughter1_pdgId)== 0){
      particles_noMotherNoDaughter->Fill(igenp.pdgId());           
      particles_noMotherNoDaughter_pT->Fill(igenp.pt());  
      particles_noMotherNoDaughter_eta->Fill(igenp.eta());  
      particles_noMotherNoDaughter_phi->Fill(igenp.phi()); 
      
    } 
  }
  
 
  NHiggs->Fill(higgs.size(),weight);
  NW->Fill(wboson.size(),weight);
  NZ->Fill(zboson.size(),weight);
 
  Nbottom->Fill(bquarks.size(),weight);
  Ntop->Fill(top.size(),weight);

  Nlept->Fill(electrons.size()+muons.size(),weight);
  Nmu->Fill(muons.size(),weight);
  Nelectrons->Fill(electrons.size());

  VLQ_Number->Fill(vlq.size(),weight);
  
  sort_by_pt(vlq);

  sort_by_pt(higgs);
  sort_by_pt(wboson);
  sort_by_pt(zboson);
  sort_by_pt(top);

  sort_by_pt(bquarks);
  sort_by_pt(electrons);
  sort_by_pt(muons);

  //vectorlike quarks
  for(auto particle: vlq){
    VLQ_pdgId->Fill(particle.pdgId());
  }
  
  if(vlq.size()>0){
    VLQ_eta_lead->Fill(vlq[0].eta(),weight);
    VLQ_phi_lead->Fill(vlq[0].phi(),weight);
    VLQ_pt_lead->Fill(vlq[0].pt(),weight);   
  }

  if(vlq.size()>1){
    VLQ_eta_subl->Fill(vlq[1].eta(),weight);
    VLQ_phi_subl->Fill(vlq[1].phi(),weight);
    VLQ_pt_subl->Fill(vlq[1].pt(),weight);   
  }
  

  //higgs plots
  for(auto particle: higgs){
    higgs_pt->Fill(particle.pt(),weight);
    higgs_eta->Fill(particle.eta(),weight);
    higgs_phi->Fill(particle.phi(),weight);
  }

  if(higgs.size()>0){
    higgs_pt_lead->Fill(higgs[0].pt(),weight); 
    higgs_eta_lead->Fill(higgs[0].eta(),weight); 
    higgs_phi_lead->Fill(higgs[0].phi(),weight); 
  }
  if(higgs.size()>1){
    higgs_pt_subl->Fill(higgs[0].pt(),weight); 
    higgs_eta_subl->Fill(higgs[0].eta(),weight); 
    higgs_phi_subl->Fill(higgs[0].phi(),weight); 
  }

  //wboson plots
  for(auto particle: wboson){
    //cout<<"pT "<<particle.pt()<<" phi "<<particle.phi()<<" eta "<< particle.eta() <<endl;
    W_pt->Fill(particle.pt(),weight);
    W_eta->Fill(particle.eta(),weight);
    W_phi->Fill(particle.phi(),weight);
  }

  if(wboson.size()>0){
    W_pt_lead->Fill(wboson[0].pt(),weight); 
    W_eta_lead->Fill(wboson[0].eta(),weight); 
    W_phi_lead->Fill(wboson[0].phi(),weight); 
  }
  if(wboson.size()>1){
    W_pt_subl->Fill(wboson[0].pt(),weight); 
    W_eta_subl->Fill(wboson[0].eta(),weight); 
    W_phi_subl->Fill(wboson[0].phi(),weight); 
  }

  //zboson plots
  for(auto particle: zboson){
    Z_pt->Fill(particle.pt(),weight);
    Z_eta->Fill(particle.eta(),weight);
    Z_phi->Fill(particle.phi(),weight);
  }

  if(zboson.size()>0){
    Z_pt_lead->Fill(zboson[0].pt(),weight); 
    Z_eta_lead->Fill(zboson[0].eta(),weight); 
    Z_phi_lead->Fill(zboson[0].phi(),weight); 
  }
  if(zboson.size()>1){
    Z_pt_subl->Fill(zboson[0].pt(),weight); 
    Z_eta_subl->Fill(zboson[0].eta(),weight); 
    Z_phi_subl->Fill(zboson[0].phi(),weight); 
  }


  //top plots
  for(auto particle: top){
    top_pt->Fill(particle.pt(),weight);
    top_eta->Fill(particle.eta(),weight);
    top_phi->Fill(particle.phi(),weight);
  }

  if(top.size()>0){
    top_pt_lead->Fill(top[0].pt(),weight); 
    top_eta_lead->Fill(top[0].eta(),weight); 
    top_phi_lead->Fill(top[0].phi(),weight); 
  }
  if(top.size()>1){
    top_pt_subl->Fill(top[0].pt(),weight); 
    top_eta_subl->Fill(top[0].eta(),weight); 
    top_phi_subl->Fill(top[0].phi(),weight); 
  }
 
 
  
}

VLQGenHists::~VLQGenHists(){}
