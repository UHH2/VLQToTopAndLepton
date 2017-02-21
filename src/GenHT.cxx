#include "UHH2/VLQToTopAndLepton/include/GenHT.h"


using namespace uhh2;
using namespace std;

GenHT::GenHT(Context & ctx){
  GenHT_val = ctx.declare_event_output<double>("GenHT");
}

bool GenHT::process(Event & event){
  if(event.isRealData) return false;
  int i=0;
  double result = 0;
  for(auto genp : *event.genparticles){
    i++;
      if(i<=2 || genp.status()<20 || genp.status()>29)continue;
      int id = abs(genp.pdgId());
      if((id >= 1 && id <= 5) || (id == 21))
	result +=genp.pt();
      //cout<<"pdg Id "<<genp.pdgId()<<" status "<<genp.status()<<endl;
    }

  event.set(GenHT_val,result);
  //cout<<result<<endl;
  return true;
}
