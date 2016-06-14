#pragma once

#include "UHH2/core/include/LorentzVector.h"

//ROOT classes
#include "TMinuit.h"
#include "TVector3.h"
#include <boost/math/constants/constants.hpp>

class NeutrinoFit{
 public:
  NeutrinoFit(bool fit = true);
  ~NeutrinoFit();
  std::vector<LorentzVector> NeutrinoFitPolar(LorentzVector lepton, LorentzVector met);
  
 private:
  static double DeltaPolarNeutrino(double PhiN, double metPx, double metPy, double PhiLep, double PtLep);
  static void polarminuitfunc(int& nDim, double* gout, double& result, double par[], int flg);
  TMinuit* positiv;
  double deltaPhiAbs(double x1, double x2);
  bool do_fit;
};
    
