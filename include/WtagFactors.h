#include "TF1.h"
#include "TFile.h"


#include "UHH2/core/include/AnalysisModule.h"


class WtagFactors: public AnalysisModule {
public:
  explicit WtagFactors(uhh2::Context & ctx, TopJetId wtag_, string jetcollname);
  virtual bool process(uhh2::Event & event) override;
private:
  TF1* puppisd_corrGEN, puppisd_corrRECO_cen, puppisd_corrRECO_for;
  uhh2::Event::Handle<double> weight;
  uhh2::Event::Handle<std::vector<TopJet>> topjet_collection;
  TopJetId wtag;
};
