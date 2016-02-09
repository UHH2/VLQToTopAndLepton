#include "UHH2/VLQToTopAndLepton/include/BprimeUncerHists.h"

using namespace std;
using namespace uhh2;

BprimeUncerHists::BprimeUncerHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyp_name): Hists(ctx, dirname){
  recohyp = ctx.get_handle<BprimeContainer>(hyp_name);

  vector<string> scale_names = {"muFup","muFdown","muRup","muRupmuFup","muRdown","muRdownmuFdown"};
  for(auto name : scale_names)
    scale.push_back(book<TH1F>("Uncertainty_scale_"+name,"B Mass[GeV]",30,50,3000));
  for(unsigned int i = 0; i<103; i++){
    pdf.push_back(book<TH1F>("Uncertainty_pdf_weight_"+to_string(i),"B Mass[GeV]",30,50,3000));
  }
}


void BprimeUncerHists::fill(const uhh2::Event & event){
  double weight = event.weight;
  BprimeContainer hyp = event.get(recohyp);
  LorentzVector whad = hyp.get_wHad();
  LorentzVector wlep = hyp.get_wLep();
  LorentzVector thad = hyp.get_topHad();
  LorentzVector tlep = hyp.get_topLep();
  int recotype = hyp.get_RecoTyp();
  LorentzVector bprime;

  if(event.isRealData) return;


  if(recotype==11||recotype==6){
    bprime = tlep+whad;
  }
  else if(recotype==12 || recotype==2){
    bprime = thad+wlep;
  }

  for(unsigned int i = 0 ; i<6; i++){
    if(event.genInfo->systweights().size()==0) break;
    unsigned int entry = scale_entries[i];
    double tmp_weight = event.genInfo->systweights().at(entry)/event.genInfo->originalXWGTUP();
    scale[i]->Fill(sqrt(bprime.M2()),weight*tmp_weight);
  }
  for(unsigned int i = 9; i<111; i++){ 
    if(event.genInfo->systweights().size()==0) break;
    double tmp_weight = event.genInfo->systweights().at(i)/event.genInfo->originalXWGTUP();
    pdf[i-9]->Fill(sqrt(bprime.M2()),weight*tmp_weight);
  }
}

BprimeUncerHists::~BprimeUncerHists(){}
