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
  vector<int> hadronicTop {6,24,-54321}; 
  vector<int> leptonicTopMu {6,24,13};
  vector<int> leptonicTopEle {6,24,11};
  vector<int> lepWMu {24,13};
  vector<int> lepWEle {24,11};
  vector<int> hadW {24,-54321};
  
  LorentzVector result(0,0,0,0);

  //top quark
  result = family(hadronicTop, event);
  myDecay.set_topHad(result);
  result = family(leptonicTopMu, event);
  if(!result.pt()>0)family(leptonicTopEle, event);
  myDecay.set_topLep(result);

  //w boson
  myDecay.set_wHad(family(hadW, event));
  result = family(lepWMu, event);
  if(!result.pt()> 0) result = family(lepWEle, event);
  myDecay.set_wLep(result);

  //VLQ && Whadd decays and b from top
  for(auto genp : *event.genparticles){
    const GenParticle* mother1 = genp.mother(event.genparticles,1);
    const GenParticle* mother2 = genp.mother(event.genparticles,2);
    int mother1_pdgId=0,mother2_pdgId=0;
    if(mother1)mother1_pdgId = mother1->pdgId();
    if(mother2)mother2_pdgId = mother2->pdgId();

    if(abs(genp.pdgId())>1000)
      myDecay.set_bprime(genp.v4());
    else if(abs(genp.pdgId())==5 && (abs(mother2_pdgId) == 6 || abs(mother1_pdgId)==6 ))
      myDecay.push_bquark(genp.v4());
    else if(abs(genp.pdgId() < 6) && (abs(mother2_pdgId) == 24 || abs(mother1_pdgId)==24))
      myDecay.push_Whadd(genp.v4());
  }

  /*
  cout<<"GenParticles ====================================="<<endl;
  cout<<"B: "<<myDecay.get_bprime().M()<<" w had: "<<myDecay.get_wHad().M()<<" w lep: "<<myDecay.get_wLep().M()<< " top lep: "<<myDecay.get_topLep().M()<<" top had: "<<myDecay.get_topHad().M()<<endl;

  cout<<"GenParticles ====================================="<<endl;  
  */
  event.set(BprimeGenLevel,move(myDecay)); 
  return true;
}


LorentzVector BprimeGen::family(vector<int> ties, const uhh2::Event & event){
  quark_decay = 0;
  for(unsigned int i =0; i< ties.size(); ++i)
    if(ties.at(i)==-54321) quark_decay = i;
  LorentzVector result = searchMother(ties,event);
  if(result.pt()==0) result = searchDaughter(ties,event);
  return result;
}

LorentzVector BprimeGen::searchMother(vector<int> familyTies ,const uhh2::Event & event){
  for(auto genp : *event.genparticles){
    if(abs(genp.pdgId()) == familyTies.at(familyTies.size()-1)){
      const GenParticle* mother1 = genp.mother(event.genparticles,1);
      const GenParticle* mother2 = genp.mother(event.genparticles,2);
      for(int m = familyTies.size()-2; m >= 0; --m){
	//cout<<familyTies.at(m)<< " search";//<<endl;
	int mother1_pdgId =0;
	int mother2_pdgId =0;
	if(mother1)mother1_pdgId = mother1->pdgId();
	if(mother2)mother2_pdgId = mother2->pdgId();
	 //cout<<" pdgId: "<< mother1_pdgId<<" "<<mother2_pdgId<<endl;
         if(abs(mother1_pdgId)!=familyTies.at(m) && abs(mother2_pdgId)!=familyTies.at(m) ){
	   bool hadronic = false;
	   if(quark_decay>0)
	     for(unsigned int i = 1; i<6; ++i){
	       if(abs(mother1_pdgId)==i){ hadronic = true; break;}
	       if(abs(mother2_pdgId)==i){ hadronic = true; break; }	
	     }
	   if(!hadronic) break;
	 }

         if(m==0){ 
	   if(abs(mother1_pdgId)==familyTies.at(0))
	     return mother1->v4();
	   else if(abs(mother2_pdgId)==familyTies.at(0))
	     return mother2->v4();
	 }
	 if(abs(mother1_pdgId)==familyTies.at(m)){
           mother1 = mother1->mother(event.genparticles,1);
           mother2 = mother1->mother(event.genparticles,2);
	 }
	 else if(abs(mother2_pdgId)==familyTies.at(m)){
	   mother1 = mother2->mother(event.genparticles,1);
           mother2 = mother2->mother(event.genparticles,2);
	 }
       }
     }
   }
  return LorentzVector(0,0,0,0);
}

LorentzVector BprimeGen::searchDaughter(vector<int> familyTies, const uhh2::Event & event){
   for(auto genp : *event.genparticles){
     if(abs(genp.pdgId()) == familyTies.at(0)){
       const GenParticle* daughter1 = genp.daughter(event.genparticles,1);
       const GenParticle* daughter2 = genp.daughter(event.genparticles,2);
       for(unsigned int m = 1; m < familyTies.size(); ++m){
         ///cout<<"search "  <<familyTies.at(m);//<<endl;
         int daughter1_pdgId =0;
         int daughter2_pdgId =0;
         if(daughter1)daughter1_pdgId = daughter1->pdgId();
         if(daughter2)daughter2_pdgId = daughter2->pdgId();
	 // cout<<" particle "<<genp.pdgId() <<" daughters: "<< daughter1_pdgId<<" "<<daughter2_pdgId<<endl;
         if(abs(daughter1_pdgId)!=familyTies.at(m) && abs(daughter2_pdgId)!=familyTies.at(m) ){
	   bool hadronic = false;
	   if(quark_decay>0)
	     for(unsigned int i = 1; i<6; ++i){
	       if(abs(daughter1_pdgId)==i){ hadronic = true; break;}
	       if(abs(daughter2_pdgId)==i){ hadronic = true; break; }	
	     }
	   if(!hadronic) break;
	 }
	 if(m==(familyTies.size()-1))
	   return genp.v4();
         if(abs(daughter1_pdgId)==familyTies.at(m)){
           daughter1 = daughter1->daughter(event.genparticles,1);
           daughter2 = daughter1->daughter(event.genparticles,2);
         }
         else if(abs(daughter2_pdgId)==familyTies.at(m)){
           daughter1 = daughter2->daughter(event.genparticles,1);
           daughter2 = daughter2->daughter(event.genparticles,2);
         }
       }
     }
   }
   //cout<<"fails"<<endl;
   return LorentzVector(0,0,0,0);
}
