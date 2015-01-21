#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"

#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"


HistFactory::HistFactory(Context& ctx):m_ctx(ctx){}

					
HistFactory::~HistFactory(){}



void HistFactory::addHists(string histClass, string histName){

  //no cut Histograms
  if(histClass.compare("ElectronHists")==0) {
    factoryHists.push_back(new ElectronHists(m_ctx,histName));
  }
  else if(histClass.compare("MuonHists")==0){
    factoryHists.push_back(new MuonHists(m_ctx,histName));
  }
  else if(histClass.compare("JetHists")==0){
    factoryHists.push_back(new JetHists(m_ctx,histName));
  }
  else if(histClass.compare("EventHists")==0){
    factoryHists.push_back(new EventHists(m_ctx,histName));
  }

  //Histograms with cuts
  for(auto cutName : cutNames){

    stringstream ss;
    ss<<cutName<<"_"<<histName;

    if(histClass.compare("ElectronHists")==0) {
      factoryHists.push_back(new ElectronHists(m_ctx,ss.str().c_str()));
    }
    else if(histClass.compare("MuonHists")==0){
      factoryHists.push_back(new MuonHists(m_ctx,ss.str().c_str()));
    }
    else if(histClass.compare("JetHists")==0){
      factoryHists.push_back(new JetHists(m_ctx,ss.str().c_str()));
    }
    else if(histClass.compare("EventHists")==0){
      factoryHists.push_back(new EventHists(m_ctx,ss.str().c_str()));
    }
  }
  
}

bool HistFactory::passAndFill(Event & event, int passOption){
  bool passCuts =true;

  for(unsigned int i = 0;i<selectionClasses.size();++i)
    factoryHists[i*selectionClasses.size()].fill(event);


  for(auto & selection : selectionClasses){
    if(selection.passes(event)){
      for(unsigned int i = 0;i<selectionClasses.size();++i)
	factoryHists[(i+1)*selectionClasses.size()].fill(event);
    }
    if(!selection.passes(event) && passOption==1) return false;
    if(!selection.passes(event)) passCuts = false;
  }

  return passCuts;

}







/*
 h_ele_jetCuts.reset(new ElectronHists(ctx, "ele_jetCuts"));
  h_jets_jetCuts.reset(new JetHists(ctx, "jets_jetCuts"));
  h_muon_jetCuts.reset(new MuonHists(ctx, "muon_jetCuts"));
  h_event_jetCuts.reset(new EventHists(ctx,"even_jetCuts"));
*/
