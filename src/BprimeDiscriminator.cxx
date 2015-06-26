#include "UHH2/VLQToTopAndLepton/include/BprimeDiscriminator.h"


using namespace std;
using namespace uhh2;


BprimeDiscriminator::BprimeDiscriminator(uhh2::Context & ctx, discriminatorType dis_, const std::string& RecoLabel,const std::string& GenLabel){
  if(RecoLabel.empty())hyps = ctx.get_handle<std::vector<BprimeContainer>>("BprimeReco");
  else  hyps = ctx.get_handle<std::vector<BprimeContainer>>(RecoLabel);
  if(GenLabel.empty())  gen = ctx.get_handle<BprimeGenContainer>("BprimeGen");
  else gen = ctx.get_handle<BprimeGenContainer>(GenLabel);
  dis=dis_;
  string disName= "DiscriminatorType_"+to_string(dis);
  //cout <<disName<<endl;
  resultHyp = ctx.declare_event_output<BprimeContainer>(disName);
}

bool BprimeDiscriminator::process(uhh2::Event & event){
  BprimeContainer result;
  if(dis ==0)
     result = ttbar_dis(event);
  else if(dis >0 && dis <4)
    result  = chiCombo_dis(event);
  else if(dis == 4)
    result  = cmsTopTag_dis(event);
  else if(dis == 5)
    result  = cmsTopTag_dis(event);

  if(result.get_RecoTyp()==-1) return false;
  //if(dis >0 && dis <4 && result.get_chiVal()>15) return false;
  //if(dis ==0 && result.get_chiVal()<44) return false;
  event.set(resultHyp,move(result));
  return true;
}

BprimeContainer BprimeDiscriminator::ttbar_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double ttbarchi =-1;
  double recoType =-1;
  for(auto hyp :  event.get(hyps)){
    LorentzVector whad = hyp.get_wHad();
    LorentzVector wlep = hyp.get_wLep();
    LorentzVector topJets = hyp.get_topJets();
    double ttbar = ((topJets+wlep).M()-174.7)*((topJets+wlep).M()-174.7)/(14*14)+(whad.M()-172.1)*(whad.M()-172.1)/(21.8*21.8);
    if(ttbar < ttbarchi || ttbarchi == -1){
      ttbarchi=ttbar;
      bestHyp =hyp;
      recoType =0;
    }
  }
  bestHyp.set_chiVal(ttbarchi);
  bestHyp.set_RecoTyp(recoType);
  if(ttbarchi ==-1) return bestHyp;
  bestHyp.set_topHad(bestHyp.get_wHad());
  bestHyp.set_topLep(bestHyp.get_wLep()+bestHyp.get_topJets());
  return bestHyp;

}

BprimeContainer BprimeDiscriminator::chiCombo_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double bprimechi =-1;
  int recoType = -1;
  for(auto hyp :  event.get(hyps)){
    LorentzVector whad = hyp.get_wHad();
    LorentzVector wlep = hyp.get_wLep();
    LorentzVector topJets = hyp.get_topJets();
    double lepTop =999999999;
    double hadTop =999999999;
    if(dis==1 || dis ==2)lepTop = ((topJets+wlep).M()-175)*((topJets+wlep).M()-175)/(14*14)+(whad.M()-82)*(whad.M()-82)/(20*20);
    if(dis==1 || dis ==3)hadTop = ((topJets+whad).M()-178)*((topJets+whad).M()-178)/(14*14)+(whad.M()-82)*(whad.M()-82)/(20*20);

    if(dis==1){
      double combochi = lepTop > hadTop ? hadTop : lepTop;
      int recoType_help = lepTop > hadTop ? 12 : 11;
      if(combochi < bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = combochi;
	recoType = recoType_help;
      }
    }
    else if(dis==2){
      if(lepTop<bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = lepTop;
	recoType = 11;
      }
    }
    else if(dis==3){
      if(hadTop<bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = hadTop;
	recoType = 12;
      }
    }
  }
  if(bprimechi ==-1) return bestHyp;
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  if(recoType==12) bestHyp.set_topHad(bestHyp.get_wHad()+bestHyp.get_topJets());
  else if(recoType==11) bestHyp.set_topLep(bestHyp.get_wLep()+bestHyp.get_topJets());
  return bestHyp;
}

BprimeContainer BprimeDiscriminator::cmsTopTag_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double bprimechi =-1;
  int recoType = -1;
  for(auto hyp :  event.get(hyps)){
    LorentzVector wlep = hyp.get_wLep();
    LorentzVector topHad = hyp.get_topHad();
    double chi = pow(topHad.M()-175,2)+ pow(deltaR(wlep,topHad)-3,2);
    if(chi<bprimechi || bprimechi==-1){
      bestHyp=hyp;
      bprimechi=chi;
      recoType=2;
    }
  }
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  return bestHyp;
}

