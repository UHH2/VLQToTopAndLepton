#include "UHH2/VLQToTopAndLepton/include/VLQControlHisto.h"



VLQControlHisto::VLQControlHisto(string name){
  HistoNames= name;
}

void VLQControlHisto::bookHistos(){
  string booking(HistoNames);
  string printing(HistoNames);

  decay    = book<TH1F>(booking.append("_decay") , printing.append(" decay modes"), 30, 0, 30);
  booking = HistoNames; printing=HistoNames;
  pt_lead  = book<TH1F>(booking.append("W_pt_lead") ,  printing.append(" p_{T}(lead) [GeV/c]"), 200, 0, 2000);
  booking = HistoNames; printing=HistoNames;
  pt_subl  = book<TH1F>(booking.append("W_pt_subl") ,  printing.append(" p_{T}(sublead) [GeV/c]"), 200, 0, 2000); 
  booking = HistoNames; printing=HistoNames;
  eta_lead = book<TH1F>(booking.append("W_eta_lead"),  printing.append(" #eta (lead)"), 40, -2.5, 2.5);
  booking = HistoNames; printing=HistoNames;
  eta_subl = book<TH1F>(booking.append("W_eta_subl"),  printing.append(" #eta (sublead)"), 40, -2.5, 2.5);
  booking = HistoNames; printing=HistoNames;
  phi_lead = book<TH1F>(booking.append("W_phi_lead"),  printing.append(" #phi (lead)"), 64, -3.2, 3.2);
  booking = HistoNames; printing=HistoNames;
  phi_subl = book<TH1F>(booking.append("W_phi_subl"),  printing.append(" #phi (sublead)"), 64, -3.2, 3.2);
  

}
  
void VLQControlHisto::FillHistos(vector<GenParticle const *> part){

  




}
