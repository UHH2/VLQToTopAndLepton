#pragma once

#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/AnalysisModule.h"


#include <boost/ptr_container/ptr_vector.hpp>
#include <string>
#include <sstream>

using namespace uhh2;
using namespace std;

class HistFactory{
 public:
  HistFactory(Context& ctx);
  ~HistFactory();

  void addSelection(Selection* selection, string cutName){selectionClasses.push_back(selection); cutNames.push_back(cutName);}

  void addHists(string histClass, string histName);

  bool passAndFill(Event & event, int passOption=0);

  void printEffiToFile(string filename){dumpFile=filename;};

 private:
  boost::ptr_vector<Selection> selectionClasses;
  boost::ptr_vector<Hists> factoryHists;

  vector<string> cutNames;
  
  vector<double> weighted_count;
  vector<int> count;
  string dumpFile;

  Context& m_ctx;

};
