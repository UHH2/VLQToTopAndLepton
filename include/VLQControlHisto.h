#pragma once


#include "UHH2/core/include/Hists.h"
#include "TH1F.h"
#include <string>




class VLQControlHisto {
 public:
  ControlHisto(string name);
  ~ControlHisto();
  
  void bookHistos();
  void FillHistos(vector<GenParticle const *> part);

 private:
  string HistoNames;
      
  TH1F* decay, pt_lead, pt_subl, eta_lead, eta_subl, phi_lead, phi_subl;

};


