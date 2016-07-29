#include "UHH2/VLQToTopAndLepton/include/BprimeGen.h"

using namespace uhh2;
using namespace std;

BprimeGen::BprimeGen(uhh2::Context & ctx, const std::string & label){
 BprimeGenLevel = ctx.declare_event_output<BprimeGenContainer>(label);
}
bool BprimeGen::process(uhh2::Event & event){
  BprimeGenContainer myDecay;
  if(event.isRealData){
    event.set(BprimeGenLevel,move(myDecay)); 
    return false;
  }
  const vector<GenParticle> * genparticles = event.genparticles;
  vector<int> hadronicTop {6,24,-54321}; 
  vector<int> leptonicTopMu {6,24,13};
  vector<int> leptonicTopEle {6,24,11};
  vector<int> lepWMu {24,13};
  vector<int> lepWEle {24,11};
  vector<int> hadW {24,-54321};
  for(auto & igenp: *genparticles){    
    if(family(hadronicTop, igenp, event))
      myDecay.set_topHad(igenp.v4());
    else if(family(leptonicTopMu, igenp, event) || family(leptonicTopEle, igenp, event))
      myDecay.set_topLep(igenp.v4());
    else if(family(hadW, igenp, event))
      myDecay.set_wHad(igenp.v4());
    else if(family(lepWMu, igenp, event)||family(lepWEle, igenp, event))
      myDecay.set_wLep(igenp.v4());
    else if(igenp.pdgId()>10000)
      myDecay.set_bprime(igenp.v4());
  }
  event.set(BprimeGenLevel,move(myDecay)); 
  return true;
}


bool BprimeGen::family(vector<int> ties, const GenParticle & part, uhh2::Event & event){

  if(abs(part.pdgId()) != ties.at(0)) return false;
  const GenParticle* daughter1 = part.daughter(event.genparticles,1);
  const GenParticle* daughter2 = part.daughter(event.genparticles,2);
  for(unsigned int m = 1; m < ties.size(); ++m){
    int daughter1_pdgId =0;
    int daughter2_pdgId =0;
    if(daughter1)daughter1_pdgId = daughter1->pdgId();
    if(daughter2)daughter2_pdgId = daughter2->pdgId();
    if(abs(daughter1_pdgId)!=ties.at(m) && abs(daughter2_pdgId)!=ties.at(m)){
      bool hadronic = false;
      if(ties.at(m)==-54321)
	for(unsigned int i = 1; i<6; ++i){
	  if(abs(daughter1_pdgId)==i){ hadronic = true; break;}
	  if(abs(daughter2_pdgId)==i){ hadronic = true; break; }	
	}
      if(!hadronic) break;
    }
    if(m==(ties.size()-1))return true;   
    if(abs(daughter1_pdgId)==ties.at(m)){
      daughter1 = daughter1->daughter(event.genparticles,1);
      daughter2 = daughter1->daughter(event.genparticles,2);
    }
    else if(abs(daughter2_pdgId)==ties.at(m)){
      daughter1 = daughter2->daughter(event.genparticles,1);
      daughter2 = daughter2->daughter(event.genparticles,2);
    }
  }
  return false;
}


