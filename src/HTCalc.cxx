#include "UHH2/VLQToTopAndLepton/include/HTCalc.h"

using namespace uhh2;
using namespace std;


namespace {
  template<typename T>
  double htCore(const vector<T> & objects,const Event & event, const boost::optional<std::function<bool (const  Jet&, const Event & )>> & object_id){
    double ht = 0.0;
    for(const auto & obj : objects){
      if(object_id){
	if((*object_id)(obj, event) )
	  ht += obj.pt();
      }
      else
	ht += obj.pt();
    }
    return ht;
  }  
}


HTCalc::HTCalc(Context & ctx, const boost::optional<JetId> & jetid, const jetTyp jetColl, const std::string & collection){
  ht_handle = ctx.get_handle<double>("HT");

  m_collection = collection;
    
  if(jetColl==jet && !collection.empty()) h_jets = ctx.get_handle<std::vector<Jet> >(collection);
  else if (jetColl==topjet && !collection.empty()) h_topjets = ctx.get_handle<std::vector<TopJet> >(collection);

  m_jetTyp = jetColl;
  m_jetid = jetid;
}

bool HTCalc::process(Event & event){
  double ht = 0.0;
  
  if( m_jetTyp==jet){
    
    vector<Jet>* jets = (!m_collection.empty()) ? &event.get(h_jets) : event.jets;
    ht = htCore<Jet>(*jets,event,m_jetid);
  }
  else if (m_jetTyp==topjet){
    vector<TopJet>* topjets = (!m_collection.empty()) ? &event.get(h_topjets) : event.topjets;
    
    ht = htCore<TopJet>(*topjets,event,m_jetid);
  }

  event.set(ht_handle, ht);

  return true;
}

/*
HTSelection::HTSelection(Context & ctx, double HTmin_): HTmin(HTmin_){
  ht = ctx.get_handle<double>("HT");
}

bool HTSelection::process(const Event & event){
  event.get_state(ht);
  return HTmin<event.get(ht);
  //return false;
}
*/
