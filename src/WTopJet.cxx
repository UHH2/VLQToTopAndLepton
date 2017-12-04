#include "UHH2/VLQToTopAndLepton/include/WTopJet.h"


bool WMass::operator()(TopJet const & topjet, uhh2::Event const & event) const {
  LorentzVector Candidate (0,0,0,0);
  for(auto & subjet : topjet.subjets()){
    Candidate += subjet.v4();
  }
  if((Candidate.M2()> min*min && Candidate.M2()<max*max) || (Candidate.M2()> min*min && max == -1) )
    return true;
  return false;
}
