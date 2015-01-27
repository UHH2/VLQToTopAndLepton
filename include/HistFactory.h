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
	      string effiFileName);
  HistFactory(Context& ctx);

  ~HistFactory();

  void addSelection(unique_ptr<Selection> selection, string cutName);
  //void addSelection(shared_ptr<Selection> selection, string cutName);


  void addHists(string histClass, string histName);

  bool passAndFill(Event & event, int passOption=0);

  void printEffiToFile(string filename){dumpFile=filename;};

 private:
  void addCounter();

  vector<unique_ptr<Selection>> selectionClasses;
  vector<unique_ptr<Hists>> factoryHists;

  vector<string> cutNames;
  
  
  vector<double> weighted_count;
  vector<double> count;
  string dumpFile;

  Context& m_ctx;
  ofstream effiFile;
  string sample;
  bool effiprint;

};




