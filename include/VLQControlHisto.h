#pragma once


#include "UHH2/core/include/Hists.h"
#include "TH1F.h"
#include <string>
#include <vector>



class VLQControlHisto {
 public:
  VLQControlHisto(string name);
  ~VLQControlHisto();
  
  void bookHistos();
  void FillHistos(vector<GenParticle const *> part);

 private:
  std::string HistoNames;
      
  TH1F* decay, pt_lead, pt_subl, eta_lead, eta_subl, phi_lead, phi_subl;

};


