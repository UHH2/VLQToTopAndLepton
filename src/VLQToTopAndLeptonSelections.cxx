#include "UHH2/VLQToTopAndLepton/include/VLQToTopAndLeptonSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>

using namespace uhh2;
using namespace std;


GenParticleFilter::GenParticleFilter(int pdgId_, int nmin_, int nmax_):pdgId(pdgId_), nmin(nmin_),nmax(nmax_)  {}




bool GenParticleFilter::passes(const uhh2::Event & event){
  
  int count = 0 ;

  for(auto genp : *event.genparticles){
    if(pdgId == abs(genp.pdgId())) count++;
    if(count>nmax) return false;
  }

  return count<=nmax && count>=nmin;
}
