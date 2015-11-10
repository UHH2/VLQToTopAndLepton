#pragma once

#include "UHH2/core/include/Selection.h"
#include "UHH2/common/include/AdditionalSelections.h"

#include "UHH2/core/include/AnalysisModule.h"

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"


//#include <boost/ptr_container/ptr_vector.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>


using namespace uhh2;
using namespace std;


class HistFactory{
 public:
  HistFactory(Context& ctx,
	      const string& effiFileName="");
  ~HistFactory();

  //HistFactory clone(string addCutname);
  void addAnalysisModule(unique_ptr<uhh2::AnalysisModule> module);
  void addSelection(unique_ptr<Selection> selection, const string& cutName);
  void addAndSelection(vector<unique_ptr<Selection>> selection, const string& cutName);
  void addOrSelection(vector<unique_ptr<Selection>> selection, const string& cutName);
  //void addSelection(shared_ptr<Selection> selection, string cutName);
  void addHists(const string& histClass, const string& histName, const std::string & hyp_name = "BprimeReco");
  void addHists(const string& histName, JetId jetid);
  void addHists(const string& histName, TopJetId topjetid);
  bool passAndFill(Event& event, int passOption=0);
  void setEffiHistName(const string& name){effiHistName=name;}

 private:
  void addCounter(); 
  void create_histos();
  unsigned int count_cuts;
  std::vector<unique_ptr<Selection>> selectionClasses;
  std::vector<string> cutNames;
  std::vector<unique_ptr<Hists>> factoryHists;
  std::vector<unique_ptr<uhh2::AnalysisModule>> AnalysisModules;
  std::vector<unsigned int> orderAnalysisModules;
  //vector<string> histNames,histClasses;
  vector<double> weighted_count;
  vector<double> count;
  Context& m_ctx;

  TH1D *cutflow_raw, *cutflow_weighted; // owned by Context
  ofstream effiFile;
  bool effiprint;
  string sample, effiHistName, effiFileName;
};




