#include "UHH2/VLQToHiggsPairProd/include/GenHists.h"
#include "UHH2/core/include/Event.h"


#include <iostream>

using namespace std;
using namespace uhh2;

namespace genhists
{
    GenParticle const * findMother (GenParticle const &, vector<GenParticle> const *);
}

using namespace genhists;

GenHists::GenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // kinematical variables 

  


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

  //top 
  top_decay    = book<TH1F>("top_decay"   , "Top decay modes", 30, 0, 30);
  top_pt_lead  = book<TH1F>("top_pt_lead" , "p_{T}^{top}(lead) [GeV/c]", 200, 0, 2000);
  top_pt_subl  = book<TH1F>("top_pt_subl" , "p_{T}^{top}(sublead) [GeV/c]", 200, 0, 2000);
  top_eta_lead = book<TH1F>("top_eta_lead", "#eta_{top}(lead)", 40, -2.5, 2.5);
  top_eta_subl = book<TH1F>("top_eta_subl", "#eta_{top}(sublead)", 40, -2.5, 2.5);
  top_phi_lead = book<TH1F>("top_phi_lead", "#phi_{top}(lead)", 64, -3.2, 3.2);
  top_phi_subl = book<TH1F>("top_phi_subl", "#phi_{top}(sublead)", 64, -3.2, 3.2);

  //Z 
  Z_decay    = book<TH1F>("Z_decay"   , "Z decay modes", 30, 0, 30);
  Z_pt_lead  = book<TH1F>("Z_pt_lead" , "p_{T}^{Z}(lead) [GeV/c]", 200, 0, 2000);
  Z_pt_subl  = book<TH1F>("Z_pt_subl" , "p_{T}^{Z}(sublead) [GeV/c]", 200, 0, 2000);
  Z_eta_lead = book<TH1F>("Z_eta_lead", "#eta_{Z}(lead)", 40, -2.5, 2.5);
  Z_eta_subl = book<TH1F>("Z_eta_subl", "#eta_{Z}(sublead)", 40, -2.5, 2.5);
  Z_phi_lead = book<TH1F>("Z_phi_lead", "#phi_{Z}(lead)", 64, -3.2, 3.2);
  Z_phi_subl = book<TH1F>("Z_phi_subl", "#phi_{Z}(sublead)", 64, -3.2, 3.2);

 
  



}


void GenHists::fill(const Event & event){
    // fill the histograms. Don't forget to always use the weight when filling:
//     double weight = event.weight;
    
    std::vector<GenParticle> const * genparticles = event.genparticles;
    
    GenParticle const * tp1 = 0;
    GenParticle const * tp2 = 0;
    GenParticle const * t1 = 0;
    GenParticle const * t2 = 0;
    GenParticle const * h1 = 0;
    GenParticle const * h2 = 0;
        
    vector<GenParticle const *> bs;
    map<double, GenParticle const *> electrons;
    map<double, GenParticle const *> muons;
    
//     int number_mu = 0;
//     int number_e = 0;
//     int number_lept = 0;
    
    for (GenParticle const & igenp : *genparticles){
        if (abs(igenp.pdgId()) == 8){
            if (!tp1) tp1 = &igenp;
            else if (igenp.pt() > tp1->pt()) {tp2 = tp1; tp1 = &igenp;}
            else if (!tp2) tp2 = &igenp;
        }
        if (abs(igenp.pdgId()) == 6){
            GenParticle const * imother = findMother(igenp, genparticles);
            if (!imother) continue;
            else if (abs(imother->pdgId()) == 8){
                if (!t1) t1 = &igenp;
                else if (igenp.pt() > t1->pt()) {t2 = t1; t1 = &igenp;}
                else if (!t2) t2 = &igenp;
            }
        }
        if (abs(igenp.pdgId()) == 25){
            GenParticle const * imother = findMother(igenp, genparticles);
            if (!imother) continue;
            else if (abs(imother->pdgId()) == 8){
                if (!h1) h1 = &igenp;
                else if (igenp.pt() > h1->pt()) {h2 = h1; h1 = &igenp;}
                else if (!h2) h2 = &igenp;
            }
        }
        if (abs(igenp.pdgId()) == 5){
            bs.push_back(&igenp);
        }
        if (abs(igenp.pdgId()) == 11){ // electron
            electrons[igenp.pt()] = &igenp;
        }
        if (abs(igenp.pdgId()) == 13){ // muon
            muons[igenp.pt()] = &igenp;
        }
    }
    
    hist("Nt")->Fill((bool)t1+(bool)t2);
    hist("NH")->Fill((bool)h1+(bool)h2);
    hist("Nb")->Fill(bs.size());
    hist("Nlept")->Fill(electrons.size()+muons.size());
    hist("Nmu")->Fill(muons.size());
    hist("Ne")->Fill(electrons.size());
    
    if (tp1) {
        hist("tpPt_lead")->Fill(tp1->pt());
        hist("eta_lead")->Fill(tp1->eta());
        hist("phi_lead")->Fill(tp1->phi());
        GenParticle const * daughter1 = tp1->daughter(genparticles);
        GenParticle const * daughter2 = tp1->daughter(genparticles, 2);
        int decay1 = 0;
        int decay2 = 0;
        if (daughter1) decay1 = abs(daughter1->pdgId());
        if (daughter2) decay2 = abs(daughter2->pdgId());
        hist("tp_decay")->Fill(decay1);
        hist("tp_decay")->Fill(decay2);
    }
    
    if (tp2) {
        hist("tpPt_subl")->Fill(tp2->pt());
        hist("eta_subl")->Fill(tp2->eta());
        hist("phi_subl")->Fill(tp2->phi());
        GenParticle const * daughter1 = tp2->daughter(genparticles);
        GenParticle const * daughter2 = tp2->daughter(genparticles, 2);
        int decay1 = 0;
        int decay2 = 0;
        if (daughter1) decay1 = abs(daughter1->pdgId());
        if (daughter2) decay2 = abs(daughter2->pdgId());
        hist("tp_decay")->Fill(decay1);
        hist("tp_decay")->Fill(decay2);
    }
    
    if (t1) {
        hist("tPt_lead")->Fill(t1->pt());
    }
    
    if (t2) {
        hist("tPt_subl")->Fill(t2->pt());
    }
    
    if (h1){
        hist("hPt_lead")->Fill(h1->pt());
        GenParticle const * daughter1 = h1->daughter(genparticles);
        GenParticle const * daughter2 = h1->daughter(genparticles, 2);
        int decay1 = 0;
        int decay2 = 0;
        if (daughter1) decay1 = abs(daughter1->pdgId());
        if (daughter2) decay2 = abs(daughter2->pdgId());
        hist("h_decay")->Fill(decay1);
        hist("h_decay")->Fill(decay2);
    }
    
    if (h2){
        hist("hPt_subl")->Fill(h2->pt());
        GenParticle const * daughter1 = h2->daughter(genparticles);
        GenParticle const * daughter2 = h2->daughter(genparticles, 2);
        int decay1 = 0;
        int decay2 = 0;
        if (daughter1) decay1 = abs(daughter1->pdgId());
        if (daughter2) decay2 = abs(daughter2->pdgId());
        hist("h_decay")->Fill(decay1);
        hist("h_decay")->Fill(decay2);
    }
    
    double bTpPt_lead = 0.;
    double bTpPt_subl = 0.;
    double bHPt_lead = 0.;
    double bHPt_subl = 0.;
    
    for (vector<GenParticle const *>::const_iterator it = bs.begin(); it != bs.end(); ++it){
        GenParticle const * ib = *it;
        GenParticle const * imother = findMother(*ib, genparticles);
        if (!imother) continue;
        if (abs(imother->pdgId()) == 8){
            if (ib->pt() > bTpPt_lead) {bTpPt_subl = bTpPt_lead; bTpPt_lead = ib->pt();}
            else if (ib->pt() > bTpPt_subl) bTpPt_subl = ib->pt();
        }
        else if (abs(imother->pdgId()) == 25){
            if (ib->pt() > bHPt_lead) {bHPt_subl = bHPt_lead; bHPt_lead = ib->pt();}
            else if (ib->pt() > bHPt_subl) bHPt_subl = ib->pt();
            for (vector<GenParticle const *>::const_iterator iit = it+1; iit != bs.end(); ++iit){
                GenParticle const * ib2 = *iit;
                GenParticle const * imother2 = findMother(*ib2, genparticles);
                if (!imother2) continue;
                if (imother->index() == imother2->index()){
                    double dR = deltaR(*ib, *ib2);
                    hist("DeltaR_bb")->Fill(dR);
                    break;
                }
            }
        }
    }
    
    if (bTpPt_lead) hist("bTpPt_lead")->Fill(bTpPt_lead);
    if (bTpPt_subl) hist("bTpPt_subl")->Fill(bTpPt_subl);
    if (bHPt_lead) hist("bHPt_lead")->Fill(bHPt_lead);
    if (bHPt_subl) hist("bHPt_subl")->Fill(bHPt_subl);
    
//     if (b1 && b2){
//         double deltaR = b1->deltaR(*b2);
// //         cout << "Test2" << endl;
//         hist("DeltaR_bb")->Fill(deltaR);
// //         cout << "Test3" << endl;
//     }
    
//     cout << "Test4" << endl;
  
  
}

GenHists::~GenHists(){}

GenParticle const * genhists::findMother (GenParticle const & igenp, vector<GenParticle> const * genparticles){
    GenParticle const * imother = igenp.mother(genparticles);
    if (!imother) imother = igenp.mother(genparticles, 2);
    return imother;
}
