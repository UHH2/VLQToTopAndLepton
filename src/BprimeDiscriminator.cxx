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
  if (dis==7)
    primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");

}

bool BprimeDiscriminator::process(uhh2::Event & event){
  BprimeContainer result;
  vector<GenParticle> tops;
  vector<GenParticle> wbosons;
  if(!event.isRealData){
    for(auto genp : *event.genparticles){
      if(fabs(genp.pdgId()==6 ))
	tops.push_back(genp);
      if(fabs(genp.pdgId()==24))
	wbosons.push_back(genp);
    }
  }
  if(dis ==0){
    result = ttbar_dis(event);
    if(!event.isRealData){
      for(auto &top : tops){
	if(deltaR(result.get_topLep(),top.v4()) < result.get_gentopdistance() || result.get_gentopdistance() ==-1)
	  result.set_gentopdistance(deltaR(result.get_topLep(),top.v4()));
	if(deltaR(result.get_topHad(),top.v4()) < result.get_genWdistance() || result.get_genWdistance() ==-1)
	  result.set_genWdistance(deltaR(result.get_topHad(),top.v4()));
      }
    }
  }
  else if((dis >0 && dis <4) || dis == 7){
    result  = chiCombo_dis(event);
  }
  else if(dis == 4 || dis == 8)
    result  = cmsTopTag_dis(event);
  else if(dis == 5)
    result  = cmsTopTag_dis(event);
  else if(dis == 6)
    result = wTag_dis(event);
  //else if(dis == 7)
  //  result = gendis(event);
  

  if(dis!=0 && !event.isRealData){
    LorentzVector leadingtop(0,0,0,0);
    LorentzVector leadingW(0,0,0,0);
    for(auto &top : tops)
      if(top.pt()> leadingtop.pt())leadingtop =top.v4();
    for(auto &wboson : wbosons)
      if(wboson.pt()> leadingW.pt())leadingW=wboson.v4();

    if(result.get_topLep().pt() > 0){
      result.set_gentopdistance(deltaR(result.get_topLep(),leadingtop));
      result.set_genWdistance(deltaR(result.get_wHad(),leadingW));
    }
    if(result.get_topHad().pt() > 0){
      result.set_gentopdistance(deltaR(result.get_topHad(),leadingtop));
      result.set_genWdistance(deltaR(result.get_wLep(),leadingW));
    }
  }

  LorentzVector forwardjet(0,0,0,0); 
  for(auto jet : *event.jets)
    if(fabs(jet.eta())>fabs(forwardjet.eta()))forwardjet = jet.v4();

  if(fabs(forwardjet.eta()) <=2.4 && dis <4){
    forwardjet = LorentzVector(0,0,0,0); 
    string unusedjets = result.get_unusedJets();
    //cout<<"Not used jets "<<unusedjets<<endl;
    for(unsigned int i= 0; i<event.jets->size(); i++){
      //cout<<!bool(unusedjets[i])<<" eta jet "<< event.jets->at(i).eta()<<endl;
      if(!unusedjets[i] && fabs(event.jets->at(i).eta())>fabs(forwardjet.eta()))
	forwardjet = event.jets->at(i).v4();
    }
  }
  else if(fabs(forwardjet.eta()) <=2.4){
    forwardjet = LorentzVector(0,0,0,0); 
    LorentzVector toptagjet = result.get_topHad();
     for(unsigned int i= 0; i<event.jets->size(); i++){
       if(deltaR(toptagjet,event.jets->at(i).v4())>0.8 && fabs(event.jets->at(i).eta())>fabs(forwardjet.eta()))
	forwardjet = event.jets->at(i).v4();
     }
  }
  result.set_forwardJet(forwardjet);
  double jetiso = -1;
  if(forwardjet.E()>0)
    for(auto jet : *event.jets)
      if((jetiso > deltaR(jet.v4(),forwardjet) || jetiso ==-1) && deltaR(jet.v4(),forwardjet)>0)jetiso = deltaR(jet.v4(),forwardjet);
  result.set_jetIso(jetiso);

  if(result.get_RecoTyp()==-1){
    if(emptyHyp){
      BprimeContainer empty_result;
      event.set(resultHyp,move(empty_result));
    }
    return false;
  }
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
    
    //this is the chi2
    double ttbar = (((topJets+wlep).M()-174.7)*((topJets+wlep).M()-174.7)/(14*14)+(whad.M()-172.1)*(whad.M()-172.1)/(21.8*21.8));
    if(ttbar < ttbarchi || ttbarchi == -1){
      ttbarchi=ttbar;
      bestHyp =hyp;
      recoType =0;
      mass = sqrt((whad+wlep+topJets).M2());
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
  BprimeGenContainer genPart;
  if(dis==7) genPart = event.get(gen);

  for(auto hyp :  event.get(hyps)){
    const LorentzVector & whad    = hyp.get_wHad();
    const LorentzVector & wlep    = hyp.get_wLep();
    const LorentzVector & topJets = hyp.get_topJets();
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

    if(topJets.pt()>0&&( dis==1 || dis ==2))lepTop = (((topJets+wlep).M()-mean_topLep)*((topJets+wlep).M()-mean_topLep)/(sigma_topLep*sigma_topLep)+(whad.M()-mean_wHad)*(whad.M()-mean_wHad)/(sigma_wHad*sigma_wHad)+pow(deltaR(whad,wlep+topJets)-mean_distance,2)/(sigma_distance*sigma_distance)+pow(whad.pt()/(wlep+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio))*0.25;
    if(dis==1 || dis ==3)hadTop = (((topJets+whad).M()-mean_topHad)*((topJets+whad).M()-mean_topHad)/(sigma_topHad*sigma_topHad)+pow(deltaR(wlep,whad+topJets)-mean_distance,2)/(sigma_distance*sigma_distance)+pow(wlep.pt()/(whad+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio))/3; 
     
    if(sqrt((wlep+topJets).M2()) < 125)
      lepTop +=1000;
    if(whad.M()>120){
      lepTop +=5000;
      hadTop +=5000; 
    }
    if(sqrt((whad+topJets).M2()) < 120)
      hadTop +=1000;
   
    if(dis==1){
      // Top lep == 11 / Top had == 12
      double combochi = lepTop > hadTop ? hadTop : lepTop;
      int recoType_help = lepTop > hadTop ? 12 : 11;
      if(combochi < bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = combochi;
	recoType = recoType_help;
	mass = sqrt((wlep+whad+topJets).M2());
      }
    }
    else if(dis==2){
      if(lepTop<bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = lepTop;
	recoType = 11;
	mass = sqrt((wlep+whad+topJets).M2());
      }
    }
    else if(dis==3){
      if(hadTop<bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = hadTop;
	recoType = 12;
	mass = sqrt((wlep+whad+topJets).M2());
      }
    }
    
    else if(dis==7){
      if(genPart.get_bprime().M()<100)continue;
      //double gen_chi = pow(genPart.get_bprime().M()-(whad+wlep+topJets).M(),2)/pow(genPart.get_bprime().M()/5,2);
      double gen_chi = 0;//pow(deltaR(genPart.get_bprime(),(whad+wlep+topJets)),2)0.25;
      //double chi2_tophad = pow((whad+topJets).M()-175,2)/64;
      //double chi2_toplep = (pow((wlep+topJets).M()-175,2)/64+pow(whad.M()-80.4,2)/100);

      //gen_chi += chi2_tophad<chi2_toplep ? chi2_tophad : chi2_toplep;
      //gen_chi = sqrt(gen_chi);
      //gen_chi /= chi2_tophad<chi2_toplep ? 2 : 3;
      
      

      if(genPart.get_wLep().pt()==0 || genPart.get_wHad().pt()==0)
	continue;
      bool toplep_gen =true;
      if ((genPart.get_topHad().pt()>0 && genPart.get_topLep().pt()>0 ) || (genPart.get_topHad().pt()==0 && genPart.get_topLep().pt()==0)){
	continue;
      }
      else if(genPart.get_topHad().pt()>0)
	toplep_gen =false;

      if(!toplep_gen){ 
	gen_chi = pow(deltaR(genPart.get_topHad(),(whad+topJets)),2)/0.16+pow((whad+topJets).M()-genPart.get_topHad().M(),2)/8;
	if(deltaR(genPart.get_topHad(),(whad+topJets))>0.4  || deltaR(genPart.get_wHad(),whad)>0.4 || deltaR(genPart.get_wLep(),wlep)>0.4|| round((genPart.get_topHad()+genPart.get_wLep()).M()) != round(genPart.get_bprime().M()))
	  gen_chi +=99999;
	//cout<<"Gen B Mass from TW "<<round((genPart.get_topHad()+genPart.get_wLep()).M())<<" from Gen particle "<<round(genPart.get_bprime().M())<<" "<< bool(round((genPart.get_topHad()+genPart.get_wLep()).M()) != round(genPart.get_bprime().M())) <<endl;      
      }
      else{ 
	gen_chi = pow(deltaR(genPart.get_topLep(),wlep+topJets),2)+pow(deltaR(genPart.get_wHad(),whad),2)+pow((wlep+topJets).M()-genPart.get_topLep().M(),2)/8+pow(whad.M()-genPart.get_wHad().M(),2)/8;
	if(deltaR(genPart.get_topLep(),(wlep+topJets))>0.4  || deltaR(genPart.get_wHad(),whad)>0.4 || deltaR(genPart.get_wLep(),wlep)>0.4 || round((genPart.get_topLep()+genPart.get_wHad()).M()) != round(genPart.get_bprime().M()))
	  gen_chi +=99999;
	//cout<<"Gen B Mass from TW "<<round((genPart.get_topLep()+genPart.get_wHad()).M())<<" from Gen particle "<<round(genPart.get_bprime().M())<<" "<<bool(round((genPart.get_topLep()+genPart.get_wHad()).M()) != round(genPart.get_bprime().M()))<<endl;
      }




      if(gen_chi<bprimechi || bprimechi ==-1){
        bestHyp = hyp;
        bprimechi = gen_chi;
	//if(chi2_tophad<chi2_toplep)recoType = 22;
	//else recoType = 21;
	if(genPart.get_topHad().pt()>0)recoType = 22;
        else if(genPart.get_topLep().pt()>0)recoType = 21;
	mass = sqrt((wlep+whad+topJets).M2());
      }
    }
       
  }
  if(bprimechi ==-1) return bestHyp;
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_Mass(mass);
  if(recoType==12 || recoType==22){
    bestHyp.set_topHad(bestHyp.get_wHad()+bestHyp.get_topJets());
  }
  else if(recoType==11 || recoType==21){
    bestHyp.set_topLep(bestHyp.get_wLep()+bestHyp.get_topJets());
  }
  bool print_bprime = false;
  if((recoType==22|| recoType==21) && bestHyp.get_topJets().pt()<=0 && print_bprime) {
    cout<<"============================================================"<<endl;
    cout<<"B mass: "<< (bestHyp.get_wLep()+bestHyp.get_wHad()+bestHyp.get_topJets()).M()<<endl;
    cout<<"W mass had: "<<bestHyp.get_wHad().M()<<" lep: "<<bestHyp.get_wLep().M()<<endl;
    cout<<"W eta had: "<<bestHyp.get_wHad().eta()<<" lep: "<<bestHyp.get_wLep().eta()<<endl;
    cout<<"top mass had: "<<(bestHyp.get_wHad()+bestHyp.get_topJets()).M()<<" lep: "<<(bestHyp.get_wLep()+bestHyp.get_topJets()).M()<<" reco: "<<bestHyp.get_RecoTyp()<<endl;
    cout<<"top eta had: "<<(bestHyp.get_wHad()+bestHyp.get_topJets()).eta()<<" lep: "<<(bestHyp.get_wLep()+bestHyp.get_topJets()).eta()<<endl;
    cout<<"chi^2: "<< bestHyp.get_chiVal()<<endl;
    cout<<"============ gen particles ============"<<endl;
    cout<<"B mass: "<<genPart.get_bprime().M()<<endl;
    cout<<"W mass had : "<<genPart.get_wHad().M()<<" lep: "<<genPart.get_wLep().M()<<endl;
    cout<<"W eta had : "<<genPart.get_wHad().eta()<<" lep: "<<genPart.get_wLep().eta()<<endl;
    cout<<"top mass had: "<<genPart.get_topHad().M()<<" lep: "<<genPart.get_topLep().M()<<endl;
    cout<<"top eta had: "<<genPart.get_topHad().eta()<<" lep: "<<genPart.get_topLep().eta()<<endl;
    cout<<"top pt had : "<<genPart.get_topHad().pt()<<" lep: "<<genPart.get_topLep().pt()<<endl;
    cout<<"W pt had : "<<genPart.get_wHad().pt()<<" lep: "<<genPart.get_wLep().pt()<<endl;
    cout<<"============================================================"<<endl;
  }
  return bestHyp;
}

BprimeContainer BprimeDiscriminator::gendis(uhh2::Event & event){
  BprimeContainer bestHyp;
  int recoType = -1;
  double mass =-1;
  BprimeGenContainer genPart = event.get(gen);

  LorentzVector whad = LorentzVector(0,0,0,0);
  LorentzVector wlep = LorentzVector(0,0,0,0);
  LorentzVector top  = LorentzVector(0,0,0,0);

  if(genPart.get_bprime().M()<100)
    return bestHyp;

  NeutrinoFit FitNeutrino;
  double wlep_min = 999999999;
  LorentzVector lep = event.get(primlep).v4();  
  for(auto  wlep_fit : FitNeutrino.NeutrinoFitPolar(lep,event.met->v4())){
    if(deltaR(genPart.get_wLep(),wlep_fit)<wlep_min){
      wlep = wlep_fit;
      wlep_min = deltaR(genPart.get_wLep(),wlep_fit);
    }
  }

  for(auto jet : *event.jets){
    if(genPart.get_wHad().pt()>0 && deltaR(genPart.get_wHad(),jet.v4())<0.4)
      whad += jet.v4();
    if(genPart.get_topLep().pt()>0 && deltaR(genPart.get_topLep(),jet.v4())<0.4)
      top += jet.v4();
    if(genPart.get_topHad().pt()>0 && deltaR(genPart.get_topHad(),jet.v4())<0.4)
      top += jet.v4();
  }
  if(genPart.get_topLep().pt()>0){
    recoType = 11;
    bestHyp.set_topLep(top+wlep);
    mass = (top+wlep+whad).M();
  }
  else if(genPart.get_topHad().pt()>0){
    recoType = 12;
    bestHyp.set_topHad(top);
    mass = (top+wlep).M();
  }

  bestHyp.set_wHad(whad);
  bestHyp.set_wLep(wlep);
  bestHyp.set_chiVal(-1);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_Mass(mass);
  
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
    
    // chi2/nodf
    double chi = (pow(topHad.M()-180,2)/23/23+ pow(deltaR(wlep,topHad)-3.1,2)/0.15/0.15)*0.5;
    if(chi<bprimechi || bprimechi==-1){
      bestHyp=hyp;
      bprimechi=chi;
      recoType=2;
      mass = sqrt((wlep+topHad).M2());
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
   
    // chi2/nodf
    double chi = (pow(wHad.M()-80,2)/23/23+ pow(deltaR(wHad,toplep)-3.1,2)/0.15/0.15)*0.5;
    if(chi<bprimechi || bprimechi==-1){
      bestHyp=hyp;
      bprimechi=chi;
      recoType=6;
      mass = sqrt((wHad+toplep).M2());
    }
  }
  if(bprimechi==-1) return bestHyp;
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_Mass(mass);
  return bestHyp;
}
