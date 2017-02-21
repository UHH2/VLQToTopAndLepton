#include "UHH2/VLQToTopAndLepton/include/WJetsReweight.h"


using namespace uhh2;
using namespace std;


WJetsReweight::WJetsReweight(Context & ctx){
  string sample = ctx.get("dataset_version");
  if(sample.find("WJet") != std::string::npos && sample.find("HT") != std::string::npos){
    work = true;
    wjets_inc = new TF1("f1","[0]*ROOT::Math::lognormal_pdf(x,[1],[2],[3])" );
    double pInc[4];pInc[0]=1.27833e+08;pInc[1]=6.64895e-02;pInc[2]=-1.34719e+00;pInc[3]=-4.83168e+01;
    wjets_ht = (TF1*)wjets_inc->Clone();
    double pHT[4];pHT[0]=1.92197e+07;pHT[1]=2.08004e+00;pHT[2]=-1.05097e+00;pHT[3]=-1.31855e+02;
    wjets_inc->SetParameters(pInc);
    wjets_ht->SetParameters(pHT);
    wpt = new TF1("linear","pol1");
    double plin[2];plin[0]=1.46292;plin[1]=-8.99954e-04;
    wpt->SetParameters(plin);

    h_wjets_reweight(ctx.declare_event_output<float>("weight_wjets_ht"));
    h_wjets_reweight_up(ctx.declare_event_output<float>("weight_wjets_ht_up"));
    h_wjets_reweight_down(ctx.declare_event_output<float>("weight_wjets_ht_down"));
  }
  else
    work = false;
}


bool WJetsReweight::process(Event & event){
  if(!work || event.isRealData) return true;
  int i=0;
  double result = 0;
  double wpt_val =0;
  for(auto genp : *event.genparticles){
    i++;
    if(abs(genp.pdgId())==24)
      wpt_val = genp.pt();
    if(i<=2 || genp.status()<20 || genp.status()>29)continue;
    int id = abs(genp.pdgId());
    if((id >= 1 && id <= 5) || (id == 21))
      result +=genp.pt();
    //cout<<"pdg Id "<<genp.pdgId()<<" status "<<genp.status()<<endl;
  }
  //cout<<"Gen HT "<<result<<" wjets inc tf1 "<<wjets_inc->Eval(result)<<" wjets ht tf1 "<<wjets_ht->Eval(result)<<" weight factor "<<wjets_inc->Eval(result)/wjets_ht->Eval(result)<<endl; 
  //event.weight *= wjets_inc->Eval(result)/wjets_ht->Eval(result);
  if(wpt_val>820){
    event.weight *= wpt->Eval(wpt_val);
    event.set(h_wjets_reweight, wpt->Eval(wpt_val))
  }
  else{
    event.weight *= 0.727;
  }
  return true;
}
