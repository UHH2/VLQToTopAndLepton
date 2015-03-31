#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/Event.h"

#include <string>

class BprimeGenContainer{
 public:
  BprimeGenContainer(){topHad.SetPxPyPzE(0,0,0,0);topLep.SetPxPyPzE(0,0,0,0);wHad.SetPxPyPzE(0,0,0,0);wLep.SetPxPyPzE(0,0,0,0);}
  LorentzVector get_wLep(){return wLep;}
  LorentzVector get_wHad(){return wHad;}
  LorentzVector get_topHad(){return topHad;}
  LorentzVector get_topLep(){return topLep;}
  void set_wLep(LorentzVector wLep_){wLep = wLep_;}
  void set_wHad(LorentzVector wHad_){wHad = wHad_;}
  void set_topHad(LorentzVector top_){topHad = top_;}
  void set_topLep(LorentzVector top_){topLep = top_;}
 private:  
  LorentzVector wLep, wHad, topHad, topLep;  
};
