#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/Event.h"

#include <string>

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
    top_jets="";
    chi =-1;
    recoTyp=-1;
  }
  
  LorentzVector get_wLep(){return wLep;}
  LorentzVector get_wHad(){return wHad;}
  LorentzVector get_topHad(){return topHad;}
  LorentzVector get_topLep(){return topLep;}
  LorentzVector get_topJets(){return topJets;}
  std::string get_wJets(){return whad_jets;}
  double get_chiVal(){return chi;}
  int get_RecoTyp(){return recoTyp;}
  
  void set_topJets(LorentzVector topJets_){topJets = topJets_;}
  void set_wLep(LorentzVector wLep_){wLep = wLep_;}
  void set_wHad(LorentzVector wHad_){wHad = wHad_;}
  void set_topLep(LorentzVector top_){topLep = top_;}
  void set_topHad(LorentzVector top_){topHad = top_;}
  void set_wJets(std::string jets){whad_jets = jets;}
  void set_chiVal(double chi_){chi=chi_;}
  void set_RecoTyp(int i_){recoTyp= i_;}
  
 private:  
  int recoTyp;
  LorentzVector wLep, wHad, topLep, topHad, topJets;  
  std::string whad_jets;
  std::string top_jets;//additional jets for the top
  Particle m_lepton;
  double chi;
};
