#include "UHH2/VLQToTopAndLepton/include/BprimeDiscriminator.h"

using namespace std;
using namespace uhh2;

BprimeDiscriminator::BprimeDiscriminator(uhh2::Context & ctx, discriminatorType dis_, const std::string& RecoLabel,const std::string Outputname, const std::string& GenLabel){
  if(RecoLabel.empty())hyps = ctx.get_handle<std::vector<BprimeContainer>>("BprimeReco");
  else  hyps = ctx.get_handle<std::vector<BprimeContainer>>(RecoLabel);
  if(GenLabel.empty())  gen = ctx.get_handle<BprimeGenContainer>("BprimeGen");
  else gen = ctx.get_handle<BprimeGenContainer>(GenLabel);
  dis=dis_;
  string disName= "DiscriminatorType_"+to_string(dis);
  emptyHyp=false;
  if(!Outputname.empty())
    disName = Outputname;
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
  else if (dis ==6)
    result = wTag_dis(event);

  if(result.get_RecoTyp()==-1){
    if(emptyHyp){
      BprimeContainer empty_result;
      event.set(resultHyp,move(empty_result));
    }
    return false;
  }
 
  /*
  else if(result.get_chiVal()> 1000){
    BprimeContainer empty_result;
    if(emptyHyp)event.set(resultHyp,move(empty_result));
    return false;
  }
  */
    
  //if(dis >0 && dis <4 && result.get_chiVal()>15) return false;
  //if(dis ==0 && result.get_chiVal()<44) return false;
  event.set(resultHyp,move(result));
  return true;
}

BprimeContainer BprimeDiscriminator::ttbar_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double ttbarchi =-1;
  double recoType =-1;
  double mass =-1;
  for(auto hyp :  event.get(hyps)){
    const LorentzVector & whad    = hyp.get_wHad();
    const LorentzVector & wlep    = hyp.get_wLep();
    const LorentzVector & topJets = hyp.get_topJets();
    mass = sqrt((whad+wlep+topJets).M2());
    //this is the chi2
    double ttbar = (((topJets+wlep).M()-174.7)*((topJets+wlep).M()-174.7)/(14*14)+(whad.M()-172.1)*(whad.M()-172.1)/(21.8*21.8));
    if(ttbar < ttbarchi || ttbarchi == -1){
      ttbarchi=ttbar;
      bestHyp =hyp;
      recoType =0;
    }
  }
  if(ttbarchi ==-1) return bestHyp;
  bestHyp.set_chiVal(ttbarchi);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_topHad(bestHyp.get_wHad());
  bestHyp.set_topLep(bestHyp.get_wLep()+bestHyp.get_topJets());
  bestHyp.set_Mass(mass);
  return bestHyp;
}

BprimeContainer BprimeDiscriminator::chiCombo_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double bprimechi =-1;
  int recoType =-1;
  double mass =-1;
  for(auto hyp :  event.get(hyps)){
    const LorentzVector & whad    = hyp.get_wHad();
    const LorentzVector & wlep    = hyp.get_wLep();
    const LorentzVector & topJets = hyp.get_topJets();

    mass = sqrt((wlep+whad+topJets).M2());
    double lepTop =999999999;
    double hadTop =999999999;
    // chi2/nodf
    //if(dis==1 || dis ==3)hadTop = (((topJets+whad).M()-178)*((topJets+whad).M()-178)/(14*14)+(whad.M()-82)*(whad.M()-82)/(14*14)+pow(deltaR(wlep,whad+topJets)-3.14,2)/0.15/0.15)*0.5;
    //if(dis==1 || dis ==2)lepTop = (((topJets+wlep).M()-175)*((topJets+wlep).M()-175)/(14*14)+(whad.M()-82)*(whad.M()-82)/(14*14)+pow(deltaR(whad,wlep+topJets)-3.14,2)/0.15/0.15)/3;
    //if(dis==1 || dis ==3)hadTop = (((topJets+whad).M()-178)*((topJets+whad).M()-178)/(14*14)+pow(deltaR(wlep,whad+topJets)-3.14,2)/0.15/0.15)*0.5;

    double mean_topLep = 175., sigma_topLep = 18.7, mean_wHad = 81., sigma_wHad = 9.6; 
    double mean_topHad = 167., sigma_topHad = 21.93;
    double mean_distance = 3.14, sigma_distance = 0.15;
    double mean_ptratio = 1, sigma_ptratio = 0.5;

    if(dis==1 || dis ==2)lepTop = (((topJets+wlep).M()-mean_topLep)*((topJets+wlep).M()-mean_topLep)/(sigma_topLep*sigma_topLep)+(whad.M()-mean_wHad)*(whad.M()-mean_wHad)/(sigma_wHad*sigma_wHad)+pow(deltaR(whad,wlep+topJets)-mean_distance,2)/(sigma_distance*sigma_distance)+pow(whad.pt()/(wlep+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio))*0.25;
    if(dis==1 || dis ==3)hadTop = (((topJets+whad).M()-mean_topHad)*((topJets+whad).M()-mean_topHad)/(sigma_topHad*sigma_topHad)+pow(deltaR(wlep,whad+topJets)-mean_distance,2)/(sigma_distance*sigma_distance)+pow(wlep.pt()/(whad+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio))/3; 
    
    if(sqrt((wlep+topJets).M2()) < 125)
      lepTop +=999999;
    if(whad.M()>120){
      lepTop +=999999;
      hadTop +=999999; 
    }

    if(sqrt((whad+topJets).M2()) < 120)
      hadTop +=1000;
    

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
  bestHyp.set_Mass(mass);
  if(recoType==12) bestHyp.set_topHad(bestHyp.get_wHad()+bestHyp.get_topJets());
  else if(recoType==11) bestHyp.set_topLep(bestHyp.get_wLep()+bestHyp.get_topJets());
  return bestHyp;
}

BprimeContainer BprimeDiscriminator::cmsTopTag_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double bprimechi =-1;
  int recoType = -1;
  double mass =-1;
  for(auto hyp :  event.get(hyps)){
    const LorentzVector & wlep = hyp.get_wLep();
    const LorentzVector & topHad = hyp.get_topHad();
    mass = sqrt((wlep+topHad).M2());
    // chi2/nodf
    double chi = (pow(topHad.M()-180,2)/23/23+ pow(deltaR(wlep,topHad)-3.1,2)/0.15/0.15)*0.5;
    if(chi<bprimechi || bprimechi==-1){
      bestHyp=hyp;
      bprimechi=chi;
      recoType=2;
    }
  }
  if(bprimechi==-1) return bestHyp;
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_Mass(mass);
  return bestHyp;
}

BprimeContainer BprimeDiscriminator::wTag_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double bprimechi =-1;
  int recoType = -1;
  double mass =-1;
  for(auto hyp :  event.get(hyps)){
    const LorentzVector & toplep = hyp.get_topLep();
    const LorentzVector & wHad = hyp.get_wHad();
    mass = sqrt((wHad+toplep).M2());
    // chi2/nodf
    double chi = (pow(wHad.M()-80,2)/23/23+ pow(deltaR(wHad,toplep)-3.1,2)/0.15/0.15)*0.5;
    if(chi<bprimechi || bprimechi==-1){
      bestHyp=hyp;
      bprimechi=chi;
      recoType=6;
    }
  }
  if(bprimechi==-1) return bestHyp;
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_Mass(mass);
  return bestHyp;
}
