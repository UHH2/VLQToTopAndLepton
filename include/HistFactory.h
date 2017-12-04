#pragma once

#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/AnalysisModule.h"

#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"

//#include <boost/ptr_container/ptr_vector.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>


using namespace uhh2;
using namespace std;


//lets be simple without templates
//since we gonna store something if it does not pass it needs to be an analysis module 
//instead of a selection :'(
class JetUncSel : public uhh2::AnalysisModule{
 public:
  explicit JetUncSel(){}

  void AddResultCollection(vector<uhh2::Event::Handle<MET>> & results){store_results=results;};
  void AddSelections(vector<unique_ptr<Selection>> selection){
    for(unsigned int i =0 ; i<selection.size();++i)
      selection_store.push_back(move(selection.at(i)));
  }
  void AddModules(vector<unique_ptr<uhh2::AnalysisModule>> module){
    for(unsigned int i =0 ; i<module.size();++i)
      module_store.push_back(move(module.at(i)));
  }

  bool process(Event & event){
    bool final = false;
    for(unsigned int i=0; i<selection_store.size();i++){
      if(selection_store[i]->passes(event)){
	final = true;
      }
      else{
	//std::cout<<"something failed, tzz met.pt =-1"<<std::endl;
	MET & tmpmet = event.get(store_results.at(i));
	tmpmet.set_pt(-1.);
	event.set(store_results[i],tmpmet);
      }     
    }
    for(unsigned int i=0; i<module_store.size();i++){
      if(module_store[i]->process(event)){
	final = true;
      }
      else{
	//std::cout<<"something failed, tzz met.pt =-1"<<std::endl;
	MET & tmpmet = event.get(store_results.at(i));
	tmpmet.set_pt(-1.);
	event.set(store_results[i],tmpmet);
      }     
    }
    //string survive = final ? "Yes" : "No";
    //cout<<"final analysis module answere: "<<survive<<endl;
    return final;
  }
  

 private:
  vector<bool> pass_vec = vector<bool>(4,false);
  vector<unique_ptr<Selection>> selection_store;
  vector<unique_ptr<uhh2::AnalysisModule>> module_store;
  vector<uhh2::Event::Handle<MET>> store_results;
};

class HistFactory{
 public:
  HistFactory(Context& ctx,
	      const string& effiFileName="");
  ~HistFactory();

  //HistFactory clone(string addCutname);
  void addAnalysisModule(unique_ptr<uhh2::AnalysisModule> module, std::string cutName );
  void addSelection(unique_ptr<Selection> selection, const string& cutName);
  void addAndSelection(vector<unique_ptr<Selection>> selection, const string& cutName);
  void addOrSelection(vector<unique_ptr<Selection>> selection, const string& cutName);
  void addAndOrSelection(vector<vector<unique_ptr<Selection>>> sel_vec, const string& cutName);
  void addJetUncSelection(vector<unique_ptr<Selection>>& selection, vector<uhh2::Event::Handle<MET>> results, const string& cutName);
  void addJetUncAnlysisModule(vector<unique_ptr<uhh2::AnalysisModule>>& module, vector<uhh2::Event::Handle<MET>> results, const string& cutName);
  void addHists(const string& histClass, const string& histName, const std::string & hyp_name = "");
  void addHists(const string& histName, JetId jetid);
  void addHists(const string& histName, TopJetId topjetid);
  void ScaleUncer();
  bool passAndFill(Event& event, int passOption=0);
  void setEffiHistName(const string& name){effiHistName=name;}

 private:
  void fillScaleUncer(uhh2::Event& event,unsigned int i);
  void addCounter(); 
  void create_histos();
  unsigned int count_cuts;
  std::vector<unique_ptr<Selection>> selectionClasses;
  std::vector<std::string> cutNames;
  std::vector<unique_ptr<Hists>> factoryHists;
  std::vector<std::vector<unique_ptr<Hists>>> factoryUncer;
  std::vector<unique_ptr<uhh2::AnalysisModule>> AnalysisModules;
  std::vector<unsigned int> orderAnalysisModules;
  std::vector<double> weighted_count;
  std::vector<double> count;
  
  uhh2::Context& m_ctx;

  std::vector<int>jetUncCount = {-1}; 
  std::vector<std::string> uncerNames;
  TH1D *cutflow_raw, *cutflow_weighted; // owned by Context
  ofstream effiFile;
  bool effiprint;
  string sample, effiHistName, effiFileName;
};




