#pragma once

#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/fwd.h"

class GenFamilySelection: public uhh2::Selection {
 public:
   explicit GenFamilySelection(std::vector<int> familyTies, int strategy=0);
   virtual bool passes(const uhh2::Event & event);
       
 private:
   bool searchMother(const uhh2::Event & event);
   bool searchDaughter(const uhh2::Event & event);
   
   std::vector<int> familyTies;
   int strategy;
   int quark_decay;


};

class GenNSelection: public uhh2::Selection {

public:
  explicit GenNSelection(int pdgId, int nmin, int nmax, int minPt, int maxPt=-1);
   virtual bool passes(const uhh2::Event & event);

private:
   int pdgId, nmin, nmax, minPt,maxPt;
};

