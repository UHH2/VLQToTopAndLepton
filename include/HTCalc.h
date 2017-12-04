#pragma once

#include <iostream>
#include <memory>


#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"


class HTCalc: public uhh2::AnalysisModule {
 public: 
  enum jetTyp {jet, topjet};
  explicit HTCalc( uhh2::Context & ctx, 
		   const boost::optional<JetId> & jetid = boost::none, 
		   const jetTyp jetColl = jet,
		   const std::string & collection = "");
  virtual bool process( uhh2::Event & event);

 private:
  std::string m_collection;
  jetTyp m_jetTyp;
  boost::optional<JetId> m_jetid;
  uhh2::Event::Handle<double> ht_handle;
  uhh2::Event::Handle<std::vector<Jet> > h_jets;
  uhh2::Event::Handle<std::vector<TopJet> > h_topjets;
};
/*
class METSelection: public uhh2::Selection {
public:
    explicit METSelection(double met_);
    virtual bool passes(const uhh2::Event & event);
private:
    double met;
};

*/

/*
class HTSelection: public uhh2::Selection{
 public:
  //enum htType{HT, HtLep};
  explicit HTSelection(uhh2::Context & ctx, double HTmin);
  virtual bool process(const uhh2::Event & event);
  
 private:
  //htType type;
  double HTmin;
  uhh2::Event::Handle<double> ht;
};
*/
