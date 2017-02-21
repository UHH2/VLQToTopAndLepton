#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/VLQToTopAndLepton/include/GenHT.h"

using namespace std;
using namespace uhh2;



class WJetsFit: public AnalysisModule {

public:
  explicit WJetsFit(Context & ctx);
  virtual bool process(Event & event) override;

private:
  std::unique_ptr<AnalysisModule> GenHT_calc;
  uhh2::Event::Handle<double> weight;
};

WJetsFit::WJetsFit(Context & ctx){
  GenHT_calc.reset(new GenHT(ctx));
  weight = ctx.declare_event_output<double>("weight");
}


bool WJetsFit::process(Event & event){
  event.set(weight,event.weight);
  return GenHT_calc->process(event);
}
// make sure the class is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(WJetsFit)
