#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"

#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"

#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"


#include "TH1D.h"



HistFactory::HistFactory(Context& ctx,
			 const string& effiFileName):m_ctx(ctx), cutflow_raw(0), cutflow_weighted(0){

  sample = ctx.get("dataset_version");
  effiFile.open(sample+string("_")+effiFileName);

  weighted_count.push_back(0);//for the count without cut
  count.push_back(0);

  effiprint = true;
 
}

HistFactory::HistFactory(Context& ctx):m_ctx(ctx), cutflow_raw(0), cutflow_weighted(0){

  effiprint = false;
}
					
HistFactory::~HistFactory(){

  if(effiprint){
    effiFile << sample << "\n" << "\n";
    //cout<<cutNames.size()<<" "<<count.size()<<endl;
    
    
    for(unsigned int i = 1; i<count.size(); i++ ){
      effiFile << cutNames[i-1] << "\n";
      effiFile << "Total\n"; 
      effiFile << "effi: "<< count[i]<<"/"<<count[0]<<" = "<<count[i]/count[0]<<" rel effi: "<<  count[i]<<"/"<<count[i-1]<<" = "<<count[i]/count[i-1]<<"\n";
      effiFile << "Weighted\n";
      effiFile << "effi: "<< weighted_count[i]<<"/"<<weighted_count[0]<<" = "<<weighted_count[i]/weighted_count[0]<<" rel effi: "<<  weighted_count[i]<<"/"<<weighted_count[i-1]<<" = "<<weighted_count[i]/weighted_count[i-1]<<"\n";
      effiFile << "\n";
    }
    
    effiFile << "---------------------------------\n";
    effiFile << "---------------------------------\n";
    effiFile.close();
  

  }
}


void HistFactory::addSelection(unique_ptr<Selection> selection, const string& cutName){
  selectionClasses.push_back(move(selection)); 
  cutNames.push_back(cutName);
  if(effiprint)addCounter();
}
/*
void HistFactory::addSelection(unique_ptr<Selection> selection, string cutName){
  selectionClasses.push_back(move(selection)); 
  cutNames.push_back(cutName);
  if(effiprint)addCounter();
}
*/

void HistFactory::addCounter(){
  
  weighted_count.push_back(0);
  count.push_back(0);
}

void HistFactory::addHists(const string& histClass, const  string& histName){
  //no cut Histograms

  unique_ptr<Hists> histTemplate;
  
  if(histClass.compare("ElectronHists")==0) {
    histTemplate.reset(new ElectronHists(m_ctx,histName.c_str()));
  }
  else if(histClass.compare("MuonHists")==0){
    histTemplate.reset(new MuonHists(m_ctx,histName.c_str()));
  }
  else if(histClass.compare("JetHists")==0){
    histTemplate.reset(new JetHists(m_ctx,histName.c_str()));
  }
  else if(histClass.compare("EventHists")==0){
    histTemplate.reset(new EventHists(m_ctx,histName.c_str()));
  }
  else if(histClass.compare("TopJetHists")==0){
    histTemplate.reset(new TopJetHists(m_ctx,histName.c_str()));
  }
  else if(histClass.compare("VLQGenHists")==0){
    histTemplate.reset(new VLQGenHists(m_ctx,histName.c_str()));
  }
  
  factoryHists.push_back(move(histTemplate));

  //Histograms with cuts
  for(auto cutName : cutNames){

    
    stringstream ss;
    ss<<cutName<<"_"<<histName;
   
    if(histClass.compare("ElectronHists")==0) {
      histTemplate.reset(new ElectronHists(m_ctx,ss.str().c_str()));
    }
    else if(histClass.compare("MuonHists")==0){
      histTemplate.reset(new MuonHists(m_ctx,ss.str().c_str()));
    }
    else if(histClass.compare("JetHists")==0){
      histTemplate.reset(new JetHists(m_ctx,ss.str().c_str()));
    }
    else if(histClass.compare("EventHists")==0){
      histTemplate.reset(new EventHists(m_ctx,ss.str().c_str()));
    }
    else if(histClass.compare("TopJetHists")==0){
      histTemplate.reset(new TopJetHists(m_ctx,ss.str().c_str()));

    }
    else if(histClass.compare("VLQGenHists")==0){
      histTemplate.reset(new VLQGenHists(m_ctx,ss.str().c_str()));
    }



    factoryHists.push_back(move(histTemplate));
   
  }
  
}

bool HistFactory::passAndFill(const Event & event, int passOption){
  bool passCuts =true;
  
  unsigned int hist_number = factoryHists.size()/(selectionClasses.size()+1);
  unsigned int cuti = 0; 

  if(!effiHistName.empty() && cutflow_raw==NULL ) create_histos();

  if(cutflow_raw){
    cutflow_raw->Fill(cuti);
    cutflow_weighted->Fill(cuti, event.weight);
  }

  if(effiprint){
    count[cuti] = count[cuti]+1;
    weighted_count[cuti] = weighted_count[cuti]+event.weight;
  }

  for(unsigned int i = 0;i<hist_number;++i)
    factoryHists[i*(selectionClasses.size()+1)]->fill(event);

  for(auto & selection : selectionClasses){
    cuti++;
    if(selection->passes(event)){
      if(cutflow_raw){
	cutflow_raw->Fill(cuti);
	cutflow_weighted->Fill(cuti, event.weight);
      }
      for(unsigned int i = 0;i<hist_number;++i)
	factoryHists[cuti+i*(selectionClasses.size()+1)]->fill(event);
      if(effiprint){
	count[cuti] = count[cuti]+1;
	weighted_count[cuti] = weighted_count[cuti]+event.weight;
      }
    }
    else{
      if(passOption==0) return false;
      else passCuts = false;
    }
  }
  return passCuts;

}

void HistFactory::create_histos(){

  cutflow_weighted = new TH1D(( effiHistName ).c_str(), ("Cutflow '" + effiHistName + "' using weights").c_str(), cutNames.size()+1, 0, cutNames.size()+1);
  cutflow_raw = new TH1D((effiHistName + "_raw").c_str(), ("Cutflow '" + effiHistName + "' unweighted").c_str(), cutNames.size()+1, 0, cutNames.size()+1);
  for(TAxis * ax : {cutflow_raw->GetXaxis(), cutflow_weighted->GetXaxis()}){
    ax->SetBinLabel(1, "all");
    for(size_t i=0; i<cutNames.size(); ++i){
      ax->SetBinLabel(i+2, cutNames.at(i).c_str());
    }
  }
  m_ctx.put("cutflow", cutflow_raw);
  m_ctx.put("cutflow", cutflow_weighted);

}
