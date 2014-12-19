#pragma once

#include "UHH2/core/include/Hists.h"

#include "TH1F.h"
/**
 *   Example class for booking and filling histograms, the new version using AnalysisModule mechanisms.
 */

class VLQGenHists: public uhh2::Hists {
 public:
  // use the same constructor arguments as Hists for forwarding:
  VLQGenHists(uhh2::Context & ctx, const std::string & dirname);
  
  virtual void fill(const uhh2::Event & ev) override;
  virtual ~VLQGenHists();
    
 private:
  TH1F *VLQ_eta_lead, *VLQ_eta_subl, *VLQ_phi_lead, *VLQ_phi_subl, *VLQ_pt_lead, *VLQ_pt_subl, *VLQ_decay;   
  TH1F *NHiggs, *NW, *NZ;
  TH1F *Nbottom, *Ntop;
  TH1F *Nlept, *Nmu, *Nelectrons;
  TH1F *higgs_decay, *DeltaR_bb, *higgs_pt_lead, *higgs_pt_subl, *higgs_eta_lead, *higgs_eta_subl, *higgs_phi_lead, *higgs_phi_subl;
  TH1F* W_decay, *W_pt_lead, *W_pt_subl, *W_eta_lead, *W_eta_subl, *W_phi_lead, *W_phi_subl;
  TH1F* top_decay, *top_pt_lead, *top_pt_subl, *top_eta_lead, *top_eta_subl, *top_phi_lead, *top_phi_subl;
  TH1F* Z_decay, *Z_pt_lead, *Z_pt_subl, *Z_eta_lead, *Z_eta_subl, *Z_phi_lead, *Z_phi_subl;




};