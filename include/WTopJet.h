#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h" 
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Selection.h"


class WMass {
public:
  explicit WMass(double min_=60,double max_=110): min(min_),max(max_){}
    
    bool operator()(const TopJet & topjet, const uhh2::Event & event) const;
    
private:
  double min, max;
};
