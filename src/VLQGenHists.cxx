#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"
#include "UHH2/core/include/Event.h"


#include <iostream>

using namespace std;
using namespace uhh2;

namespace genhists
{
    GenParticle const * findMother (GenParticle const &, vector<GenParticle> const *);
}

using namespace genhists;

VLQGenHists::VLQGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // kinematical variables 

  
  VLQ_Number =  book<TH1F>("VLQ_number", "Number of VLQ", 5, -0.5, 4.5);
  VLQ_pdgId =  book<TH1F>("VLQ_pdgId", "pdgId", 5, 4.5, 9.5); 

  VLQ_eta_lead = book<TH1F>("VLQ_eta_lead", "#eta_{VLQ}(lead)", 40, -2.5, 2.5);
  VLQ_eta_subl = book<TH1F>("VLQ_eta_subl", "#eta_{VLQ}(sublead)", 40, -2.5, 2.5);
  
  VLQ_phi_lead = book<TH1F>("VLQ_phi_lead", "#phi_{VLQ}(lead)", 64, -3.2, 3.2);
  VLQ_phi_subl = book<TH1F>("VLQ_phi_subl", "#phi_{VLQ}(sublead)", 64, -3.2, 3.2);
  
  VLQ_pt_lead  =  book<TH1F>("VLQ_pt_lead" , "p_{T}^{VLQ}(lead) [GeV/c]", 200, 0, 2000);
  VLQ_pt_subl  =  book<TH1F>("VLQ_pt_subl" , "p_{T}^{VLQ}(sublead) [GeV/c]", 200, 0, 2000);
  
  VLQ_decay    = book<TH1F>("VLQ_decay", "VLQ decay modes", 30, 0, 30);


  //number of bosons 
  NHiggs = book<TH1F>("NHiggs" , "Number of Higgses", 5, 0, 5);
  NW     = book<TH1F>("NW"     , "Number of Higgses", 5, 0, 5);
  NZ     = book<TH1F>("NZ"     , "Number of Higgses", 5, 0, 5);

  //number of quarks
  Nbottom = book<TH1F>("Nbottom", "Number of bs", 8, 0, 8);
  Ntop = book<TH1F>("Ntop"   , "Number of tops", 5, 0, 5);

  //number of leptons
  Nlept      = book<TH1F>("Nlept"     , "Number of leptons", 15, 0, 15);
  Nmu        = book<TH1F>("Nmu"       , "Number of muons", 15, 0, 15);
  Nelectrons = book<TH1F>("Nelectrons", "Number of electrons", 15, 0, 15);
  
  //higgs 
  higgs_decay    = book<TH1F>("higgs_decay"     , "Higgs decay modes", 30, 0, 30);
  DeltaR_bb      = book<TH1F>("DeltaR_bb"   , "#Delta R_{bb}", 50, 0, 5); //b variable
  higgs_pt_lead  = book<TH1F>("higgs_pt_lead" , "p_{T}^{higgs}(lead) [GeV/c]", 200, 0, 2000);
  higgs_pt_subl  = book<TH1F>("higgs_pt_subl" , "p_{T}^{higgs}(sublead) [GeV/c]", 200, 0, 2000);
  higgs_eta_lead = book<TH1F>("higgs_eta_lead", "#eta_{higgs}(lead)", 40, -2.5, 2.5);
  higgs_eta_subl = book<TH1F>("higgs_eta_subl", "#eta_{higgs}(sublead)", 40, -2.5, 2.5);
  higgs_phi_lead = book<TH1F>("higgs_phi_lead", "#phi_{higgs}(lead)", 64, -3.2, 3.2);
  higgs_phi_subl = book<TH1F>("higgs_phi_subl", "#phi_{higgs}(sublead)", 64, -3.2, 3.2);

  //W
  W_decay    = book<TH1F>("W_decay"   , "W decay modes", 30, 0, 30);
  W_pt_lead  = book<TH1F>("W_pt_lead" , "p_{T}^{W}(lead) [GeV/c]", 200, 0, 2000);
  W_pt_subl  = book<TH1F>("W_pt_subl" , "p_{T}^{W}(sublead) [GeV/c]", 200, 0, 2000);
  W_eta_lead = book<TH1F>("W_eta_lead", "#eta_{W}(lead)", 40, -2.5, 2.5);
  W_eta_subl = book<TH1F>("W_eta_subl", "#eta_{W}(sublead)", 40, -2.5, 2.5);
  W_phi_lead = book<TH1F>("W_phi_lead", "#phi_{W}(lead)", 64, -3.2, 3.2);
  W_phi_subl = book<TH1F>("W_phi_subl", "#phi_{W}(sublead)", 64, -3.2, 3.2);

  //Z 
  Z_decay    = book<TH1F>("Z_decay"   , "Z decay modes", 30, 0, 30);
  Z_pt_lead  = book<TH1F>("Z_pt_lead" , "p_{T}^{Z}(lead) [GeV/c]", 200, 0, 2000);
  Z_pt_subl  = book<TH1F>("Z_pt_subl" , "p_{T}^{Z}(sublead) [GeV/c]", 200, 0, 2000);
  Z_eta_lead = book<TH1F>("Z_eta_lead", "#eta_{Z}(lead)", 40, -2.5, 2.5);
  Z_eta_subl = book<TH1F>("Z_eta_subl", "#eta_{Z}(sublead)", 40, -2.5, 2.5);
  Z_phi_lead = book<TH1F>("Z_phi_lead", "#phi_{Z}(lead)", 64, -3.2, 3.2);
  Z_phi_subl = book<TH1F>("Z_phi_subl", "#phi_{Z}(sublead)", 64, -3.2, 3.2);

  //top 
  top_decay    = book<TH1F>("top_decay"   , "Top decay modes", 30, 0, 30);
  top_pt_lead  = book<TH1F>("top_pt_lead" , "p_{T}^{top}(lead) [GeV/c]", 200, 0, 2000);
  top_pt_subl  = book<TH1F>("top_pt_subl" , "p_{T}^{top}(sublead) [GeV/c]", 200, 0, 2000);
  top_eta_lead = book<TH1F>("top_eta_lead", "#eta_{top}(lead)", 40, -2.5, 2.5);
  top_eta_subl = book<TH1F>("top_eta_subl", "#eta_{top}(sublead)", 40, -2.5, 2.5);
  top_phi_lead = book<TH1F>("top_phi_lead", "#phi_{top}(lead)", 64, -3.2, 3.2);
  top_phi_subl = book<TH1F>("top_phi_subl", "#phi_{top}(sublead)", 64, -3.2, 3.2);
  



}


void VLQGenHists::fill(const Event & event){
  // fill the histograms. Don't forget to always use the weight when filling:
  //     double weight = event.weight;
  double weight = event.weight;

   
  std::vector<GenParticle> const * genparticles = event.genparticles;
  
  vector<GenParticle const *> vlq;

  vector<GenParticle const *> higgs;
  vector<GenParticle const *> wboson;
  vector<GenParticle const *> zboson;

  vector<GenParticle const *> top;
  vector<GenParticle const *> bquarks;
  vector<GenParticle const *> electrons;
  vector<GenParticle const *> muons;
  


  for(auto igenp : *genparticles){
    if (abs(igenp.pdgId()) == 8 || abs(igenp.pdgId()) == 7)
      vlq.push_back(&igenp);
    else if(abs(igenp.pdgId()) == 5 )
      bquarks.push_back(&igenp);
    else if(abs(igenp.pdgId()) == 6 )
      top.push_back(&igenp);
    else if (abs(igenp.pdgId()) == 11 )
      electrons.push_back(&igenp);
    else if (abs(igenp.pdgId()) == 13 )
      muons.push_back(&igenp);
    else if(abs(igenp.pdgId()) == 23 )
      zboson.push_back(&igenp);
    else if(abs(igenp.pdgId()) == 24 )
      wboson.push_back(&igenp);
    else if(abs(igenp.pdgId()) == 25 )
      higgs.push_back(&igenp);
  }


  NHiggs->Fill(higgs.size(),weight);
  NW->Fill(wboson.size(),weight);
  NZ->Fill(zboson.size(),weight);

  Nbottom->Fill(bquarks.size(),weight);
  Ntop->Fill(top.size(),weight);

  Nlept->Fill(electrons.size()+muons.size(),weight);
  Nmu->Fill(muons.size(),weight);


  VLQ_Number->Fill(vlq.size(),weight);
  
  for(auto particle: vlq){
    VLQ_pdgId->Fill(particle.pdgId());
    
  }
  
}

VLQGenHists::~VLQGenHists(){}

GenParticle const * genhists::findMother (GenParticle const & igenp, vector<GenParticle> const * genparticles){
    GenParticle const * imother = igenp.mother(genparticles);
    if (!imother) imother = igenp.mother(genparticles, 2);
    return imother;
}
