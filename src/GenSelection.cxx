
#include "UHH2/VLQToTopAndLepton/include/GenSelection.h"


using namespace uhh2;
using namespace std;

GenFamilySelection::GenFamilySelection(std::vector<int> familyTies_, int strategy_):familyTies(familyTies_), strategy(strategy_){} 


bool GenFamilySelection::passes(const uhh2::Event & event){

    if(strategy==0)
      return searchMother(event) && searchDaughter(event);
    else if(strategy==1)
      return searchMother(event);
    else if(strategy==2)
      return searchDaughter(event);
    else if(strategy==3)
      return searchMother(event) || searchDaughter(event);
    else{
      cout<<"No strategy was choosen, aborting run"<<endl;
      assert(0==1);
    }

    return false;
}

bool GenFamilySelection::searchMother(const uhh2::Event & event){

   for(auto genp : *event.genparticles){
     if(genp.pdgId() == familyTies.at(familyTies.size()-1)){
       const GenParticle* mother1 = genp.mother(event.genparticles,1);
       const GenParticle* mother2 = genp.mother(event.genparticles,2);
       for(int m = familyTies.size()-2; m >= 0; --m){
	//cout<<familyTies.at(m)<< " search";//<<endl;
	 int mother1_pdgId =0;
         int mother2_pdgId =0;
         if(mother1)mother1_pdgId = mother1->pdgId();
	 if(mother2)mother2_pdgId = mother2->pdgId();
	 //cout<<" pdgId: "<< mother1_pdgId<<" "<<mother2_pdgId<<endl;

         if(abs(mother1_pdgId)!=familyTies.at(m) && abs(mother2_pdgId)!=familyTies.at(m) ) break;
         if(m==0) return true;
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
   return false;
}

bool GenFamilySelection::searchDaughter(const uhh2::Event & event){

   for(auto genp : *event.genparticles){
     if(genp.pdgId() == familyTies.at(0)){
       const GenParticle* daughter1 = genp.daughter(event.genparticles,1);
       const GenParticle* daughter2 = genp.daughter(event.genparticles,2);
       for(int m = 1; m < familyTies.size(); ++m){
         cout<<familyTies.at(m)<< " search";//<<endl;
         int daughter1_pdgId =0;
         int daughter2_pdgId =0;
         if(daughter1)daughter1_pdgId = daughter1->pdgId();
         if(daughter2)daughter2_pdgId = daughter2->pdgId();
         cout<<" pdgId: "<< daughter1_pdgId<<" "<<daughter2_pdgId<<endl;

         if(abs(daughter1_pdgId)!=familyTies.at(m) && abs(daughter2_pdgId)!=familyTies.at(m) ) break;
         if(m==(familyTies.size()-1)) return true;
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
   return false;
}



