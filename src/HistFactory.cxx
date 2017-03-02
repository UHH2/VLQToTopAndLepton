#include "UHH2/VLQToTopAndLepton/include/HistFactory.h"

#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/GenJetsHists.h"

#include "UHH2/VLQToTopAndLepton/include/BprimeRecoHists.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeHypHists.h"
#include "UHH2/VLQToTopAndLepton/include/BprimeUncerHists.h"
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"
#include "UHH2/VLQToTopAndLepton/include/EventKinematicHists.h"

#include "TH1D.h"
#include <math.h> 


HistFactory::HistFactory(Context& ctx,
			 const string& effiFileName_):m_ctx(ctx), cutflow_raw(0), cutflow_weighted(0){

  sample = ctx.get("dataset_version");
  effiFileName=effiFileName_;

  cutNames.push_back("");
  count_cuts =0;
  if(!effiFileName.empty()){
    effiFile.open(sample+string("_")+effiFileName);
    effiprint = true;
    weighted_count.push_back(0);//for the count without cut
    count.push_back(0);
  }
  else effiprint = false;
  
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
  count_cuts++;
}
void HistFactory::addAndSelection(vector<unique_ptr<Selection>> selection, const string& cutName){
  unique_ptr<AndSelection> myAndSel;
  myAndSel.reset(new AndSelection(m_ctx));
  for(unsigned int i=0; i<selection.size(); i++)
    myAndSel->add(to_string(i),move(selection.at(i)));
  addSelection(move(myAndSel),cutName);
}
void HistFactory::addOrSelection(vector<unique_ptr<Selection>> selection, const string& cutName){
  unique_ptr<OrSelection> myOrSel;
  myOrSel.reset(new OrSelection());
  for(unsigned int i=0; i<selection.size(); i++)
    myOrSel->add(move(selection.at(i)));
  addSelection(move(myOrSel),cutName);
}
void HistFactory::addCounter(){
  weighted_count.push_back(0);
  count.push_back(0);
}

void HistFactory::addAnalysisModule(unique_ptr<uhh2::AnalysisModule> module){
  AnalysisModules.push_back(move(module));
  orderAnalysisModules.push_back(selectionClasses.size());
}

void HistFactory::addHists(const string& histClass, const string& histName, const string & hyp_name){
  unique_ptr<Hists> histTemplate;
  vector<unique_ptr<Hists>> uncerHistsTemplate;
  for(const auto & cutName : cutNames){
    stringstream ss;
    if(!cutName.empty()) ss<<cutName<<"_"<<histName;
    else  ss<<histName;
    if(histClass.compare("ElectronHists")==0) {
      histTemplate.reset(new ElectronHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new ElectronHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("MuonHists")==0){
      histTemplate.reset(new MuonHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new MuonHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("JetHists")==0){
      histTemplate.reset(new JetHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new JetHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("GenJetHists")==0){
      histTemplate.reset(new GenJetsHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new GenJetsHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("EventHists")==0){
      histTemplate.reset(new EventHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new EventHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("TopJetHists")==0){
      histTemplate.reset(new TopJetHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new TopJetHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("VLQGenHists")==0){
      histTemplate.reset(new VLQGenHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new VLQGenHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("BprimeRecoHists")==0){
      histTemplate.reset(new BprimeRecoHists(m_ctx,ss.str().c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new BprimeRecoHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("BprimeHypHists")==0){
      histTemplate.reset(new BprimeHypHists(m_ctx,ss.str().c_str(),hyp_name));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new BprimeHypHists(m_ctx,(ss.str()+"_"+name).c_str(),hyp_name));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("BprimeUncerHists")==0){
      histTemplate.reset(new BprimeUncerHists(m_ctx,ss.str().c_str(),hyp_name));
    }
    else if(histClass.compare("LuminosityHists")==0){
      histTemplate.reset(new LuminosityHists(m_ctx,(ss.str()+"_perlumibin").c_str()));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new LuminosityHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else if(histClass.compare("EventKinematicHists")==0){
      histTemplate.reset(new EventKinematicHists(m_ctx,ss.str().c_str(),hyp_name));
      for(auto & name : uncerNames){
	unique_ptr<Hists> uncerHists;
	uncerHists.reset(new EventKinematicHists(m_ctx,(ss.str()+"_"+name).c_str(),hyp_name));
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    }
    else
      cerr<<"You ask for a not supported hist class, please check spelling or add class";
    factoryHists.push_back(move(histTemplate));
    factoryUncer.push_back(move(uncerHistsTemplate));
  }
}
void HistFactory::addHists(const string& histName, JetId jetid){
  unique_ptr<JetHists> histTemplate;
  vector<unique_ptr<Hists>> uncerHistsTemplate;
  for(const auto & cutName : cutNames){
    stringstream ss;
    if(!cutName.empty()) ss<<cutName<<"_"<<histName;
    else  ss<<histName;
    histTemplate.reset(new JetHists(m_ctx,ss.str().c_str()));
    histTemplate->set_JetId(jetid);
    factoryHists.push_back(move(histTemplate));
    for(auto & name : uncerNames){
	unique_ptr<JetHists> uncerHists;
	uncerHists.reset(new JetHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHists->set_JetId(jetid);
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    factoryUncer.push_back(move(uncerHistsTemplate));
  }
}

void HistFactory::addHists(const string& histName, TopJetId topjetid){
  unique_ptr<TopJetHists> histTemplate;
  vector<unique_ptr<Hists>> uncerHistsTemplate;
  for(const auto & cutName : cutNames){
    stringstream ss;
    if(!cutName.empty()) ss<<cutName<<"_"<<histName;
    else  ss<<histName;
    histTemplate.reset(new TopJetHists(m_ctx,ss.str().c_str()));
    histTemplate->set_TopJetId(topjetid);
    factoryHists.push_back(move(histTemplate));
    for(auto & name : uncerNames){
	unique_ptr<TopJetHists> uncerHists;
	uncerHists.reset(new TopJetHists(m_ctx,(ss.str()+"_"+name).c_str()));
	uncerHists->set_TopJetId(topjetid);
	uncerHistsTemplate.push_back(move(uncerHists));
      }
    factoryUncer.push_back(move(uncerHistsTemplate));
  }
}

//passOption: 0 event has to pass all cuts, 1 event passes this cut
bool HistFactory::passAndFill(Event & event, int passOption){
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
    for(auto & module : orderAnalysisModules)
      if(module == cuti){
	  bool modulePass = AnalysisModules.at(&module-&orderAnalysisModules[0])->process(event);
	  if(!modulePass && passOption==0) return false;
	  else if(!modulePass && passOption==1) passCuts = false;  
      }	
    cuti++;
    if(selection->passes(event)){
      if(cutflow_raw){
	cutflow_raw->Fill(cuti);
	cutflow_weighted->Fill(cuti, event.weight);
      }
      for(unsigned int i = 0;i<hist_number;++i){
	factoryHists[cuti+i*(selectionClasses.size()+1)]->fill(event);
	fillScaleUncer(event,cuti+i*(selectionClasses.size()+1));
      }
      if(effiprint){
	count[cuti] = count[cuti]+1;
	weighted_count[cuti] = weighted_count[cuti]+event.weight;
      }
    }
    else{
      if(passOption==0) return false;
      else if(passOption==1) passCuts = false;
    }
  }
  for(auto & module : orderAnalysisModules)
    if(module == cuti)  AnalysisModules[&module-&orderAnalysisModules[0]]->process(event);
  return passCuts;
}

void HistFactory::create_histos(){
  cutflow_weighted = new TH1D(( effiHistName ).c_str(), ("Cutflow '" + effiHistName + "' using weights").c_str(), cutNames.size(), 0, cutNames.size());
  cutflow_raw = new TH1D((effiHistName + "_raw").c_str(), ("Cutflow '" + effiHistName + "' unweighted").c_str(), cutNames.size(), 0, cutNames.size());
  for(TAxis * ax : {cutflow_raw->GetXaxis(), cutflow_weighted->GetXaxis()}){
    ax->SetBinLabel(1, "all");
    for(size_t i=1; i<cutNames.size(); ++i){
      ax->SetBinLabel(i+1, cutNames.at(i).c_str());
    }
  }
  m_ctx.put("cutflow", cutflow_raw);
  m_ctx.put("cutflow", cutflow_weighted);
}

//not used since histograms are needed maybe later on an optiont
void HistFactory::fillScaleUncer(Event & event,unsigned int i){
  if(event.isRealData) return;
  if(event.genInfo->systweights().size() ==0) return;
  if(uncerNames.size()==0) return;
  double syst_downWeight =1.;
  double syst_upWeight =1.;
  double old_weight = event.weight;
  for(unsigned int it=1; it <9; it++){
    if(it == 5 || it == 7) continue;
    double temp_weight = event.genInfo->systweights().at(it)/event.genInfo->originalXWGTUP();
    if(temp_weight<syst_downWeight) syst_downWeight = temp_weight;
    if(temp_weight>syst_upWeight) syst_upWeight = temp_weight;
  }
  //for(unsigned int i= 0; i<uncerNames.size();i++ ){
  //cout<<"Up weight "<<syst_upWeight<<" Down weight "<<syst_downWeight<<endl;

  event.weight *=syst_upWeight;
  factoryUncer.at(i).at(0)->fill(event);
  event.weight = old_weight;
  event.weight *= syst_downWeight;
  factoryUncer.at(i).at(1)->fill(event);
  event.weight = old_weight;
  //}
}

void HistFactory::ScaleUncer(){ 
  uncerNames.push_back("ScaleUp_plus");
  uncerNames.push_back("ScaleDown_minus");
}
