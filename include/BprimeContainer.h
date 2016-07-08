#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/Event.h"

#include <string>
#include <vector>

class BprimeContainer{
 public:
  BprimeContainer(){
    LorentzVector emptyLorentz(0,0,0,0);
    topHad = emptyLorentz;
    topLep = emptyLorentz;
    wHad = emptyLorentz;
    wLep = emptyLorentz;
    topJets = emptyLorentz;
    whad_jets="";
    unused_jets="";
    chi =-1;
    mass =-1;
    recoTyp=-1;
    btag_discriminand=-1;
    btag_whad_discr =-1;
  }
  
  LorentzVector get_wLep(){return wLep;}
  LorentzVector get_wHad(){return wHad;}
  std::vector<LorentzVector> get_wHadLorentz(){return wHad_lorentz;}
  LorentzVector get_topHad(){return topHad;}
  LorentzVector get_topLep(){return topLep;}
  LorentzVector get_topJets(){return topJets;}
  std::vector<LorentzVector> get_topLorentz(){return top_lorentz;}
  std::string get_wJets(){return whad_jets;}
  std::string get_unusedJets(){return unused_jets;}
  double get_chiVal(){return chi;}
  double get_Mass(){return mass;}
  int get_RecoTyp(){return recoTyp;}
  double get_btag_discriminator(){return btag_discriminand;}  
  double get_btag_whad(){return btag_whad_discr;}
  
  void set_btag_whad(double b_whad){btag_whad_discr =b_whad;}
  void set_btag_discriminator(double btag){btag_discriminand = btag;}
  void set_topJets(LorentzVector topJets_){topJets = topJets_;}
  void set_topLorentz(std::vector<LorentzVector> top_LorentzJets){top_lorentz = top_LorentzJets;}
  void set_wLep(LorentzVector wLep_){wLep = wLep_;}
  void set_wHad(LorentzVector wHad_){wHad = wHad_;}
  void set_wHadJets(std::vector<LorentzVector> whadJets){wHad_lorentz = whadJets;}
  void set_topLep(LorentzVector top_){topLep = top_;}
  void set_topHad(LorentzVector top_){topHad = top_;}
  //void set_forwardJet(LorentzVector forwardJet_){forwardJet = forwardJet_;}
  //void set_balanceJet(LorentzVector balanceJet_){balanceJet = balanceJet_;}
  void set_wJets(std::string jets){whad_jets = jets;}
  void set_unusedJets(std::string jets){unused_jets = jets;}
  void set_chiVal(double chi_){chi=chi_;}
  void set_Mass(double mass_){mass=mass_;}
  void set_RecoTyp(int i_){recoTyp= i_;}
  
 private:  
  int recoTyp;
  LorentzVector wLep, wHad, topLep, topHad, topJets;//, forwardJet, balanceJet;  
  std::vector<LorentzVector> wHad_lorentz, top_lorentz;
  std::string whad_jets;
  std::string unused_jets;//additional jets for the top
  Particle m_lepton;
  double chi,mass;
  double btag_discriminand, btag_whad_discr;
};
