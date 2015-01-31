#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/ObjectIdUtils.h"


class GenParticleFilter: public uhh2::Selection {
 public:
  explicit GenParticleFilter(int pdgId, int nmin, int nmax=-1);
  virtual bool passes(const uhh2::Event & event);
 private:
  int pdgId, nmin, nmax;
  

};


//class Bprime
