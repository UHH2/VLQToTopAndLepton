#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/Event.h"

#include <string>

class BprimeGenContainer{
 public:
  BprimeGenContainer(){
    topHad.SetPxPyPzE(0,0,0,0);
    topLep.SetPxPyPzE(0,0,0,0);
    wHad.SetPxPyPzE(0,0,0,0);
    wLep.SetPxPyPzE(0,0,0,0);
    bprime.SetPxPyPzE(0,0,0,0);
    associated.SetPxPyPzE(0,0,0,0);
    lepton.SetPxPyPzE(0,0,0,0);
    forward.SetPxPyPzE(0,0,0,0);
    bquark = std::vector<LorentzVector>();
    whad_daughters = std::vector<LorentzVector>();
  }
  LorentzVector get_wLep(){return wLep;}
  LorentzVector get_wHad(){return wHad;}
  LorentzVector get_topHad(){return topHad;}
  LorentzVector get_topLep(){return topLep;}
  LorentzVector get_bprime(){return bprime;}
  
  void set_wLep(LorentzVector wLep_){wLep = wLep_;}
  void set_wHad(LorentzVector wHad_){wHad = wHad_;}
  void set_topHad(LorentzVector top_){topHad = top_;}
  void set_topLep(LorentzVector top_){topLep = top_;}
  void set_bprime(LorentzVector bprime_){bprime=bprime_;}
  void set_associated(LorentzVector ass_){associated=ass_;}
  void set_forward(LorentzVector for_){forward=for_;}
  void push_bquark(LorentzVector b_){bquark.push_back(b_);}
  void push_Wlepd(LorentzVector d_){wlep_daughters.push_back(d_);}
  void push_Whadd(LorentzVector d_){whad_daughters.push_back(d_);}
  
 private:  
  LorentzVector wLep, wHad, topHad, topLep, bprime, forward, associated,  lepton;
  std::vector<LorentzVector> wlep_daughters, whad_daughters,bquark;
};
