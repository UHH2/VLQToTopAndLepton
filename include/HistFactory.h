#pragma once

#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/AnalysisModule.h"

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

  void addSelection(unique_ptr<Selection> selection, const string& cutName);
  //void addSelection(shared_ptr<Selection> selection, string cutName);
  void addHists(const string& histClass, const string& histName);
  bool passAndFill(const Event& event, int passOption=0);
  void setEffiHistName(const string& name){effiHistName=name;}

 private:
  void addCounter(); 
  void create_histos();

  vector<unique_ptr<Selection>> selectionClasses;
  vector<string> cutNames;
  vector<unique_ptr<Hists>> factoryHists;
  //vector<string> histNames,histClasses;
  vector<double> weighted_count;
  vector<double> count;
  Context& m_ctx;

  TH1D * cutflow_raw, * cutflow_weighted; // owned by Context
  ofstream effiFile;
  bool effiprint;
  string sample, effiHistName, effiFileName;
};




