#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/Event.h"

#include <string>

class BprimeContainer{
 public:
  /*
  BprimeContainer(const BprimeContainer & old){
    wLep = old.get_wLep();
    wHad = old.get_wHad();
    wHad_jets = old.get_wJets();
    }*/
  //~BprimeContainer(){};
  
  LorentzVector get_wLep(){return wLep;}
  LorentzVector get_wHad(){return wHad;}
  LorentzVector get_top(){return top;}
  LorentzVector get_topJets(){return topJets;}
  std::string get_wJets(){return whad_jets;}
  //std::vector<Particle> get_wHadJets(){return wHad_jets;}
  //std::vector<Particle> get_wLepJets(){return wLep_jets;}
  //std::vector<Particle> get_bJets(){return b_jets;}
  
  void set_topJets(LorentzVector topJets_){topJets = topJets_;}
  void set_wLep(LorentzVector wLep_){wLep = wLep_;}
  void set_wHad(LorentzVector wHad_){wHad = wHad_;}
  void set_top(LorentzVector top_){top = top_;}
  void set_wJets(std::string jets){whad_jets = jets;}
  //void set_topJets(std::string jets){top_jets = jets;}
  //void add_wLepJet(Particle wJet){wLep_jets.push_back(wJet);}
  //void add_wHadJet(Particle wJet){wHad_jets.push_back(wJet);}
  //void add_top_candiJets(Jet jetCandi){top_candi.emplace_back(jetCandi);}
  
 private:  
  LorentzVector wLep, wHad, top, topJets;  
  //b_jets refers to the b-Jet of the top decay!
  //std::vector<Jet> top_candi;
  std::string whad_jets;
  std::string top_jets;//additional jets for the top
  Particle m_lepton;
};
