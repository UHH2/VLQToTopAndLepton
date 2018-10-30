#include "UHH2/VLQToTopAndLepton/include/BprimeDiscriminator.h"


#include <boost/math/distributions/chi_squared.hpp>


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
  int forwardnumber =0;
  int forwardnumber_eta2 =0;
  
  //for(auto jet : *event.jets)
  //  if(fabs(jet.eta())>fabs(forwardjet.eta()))forwardjet = jet.v4();

  if(dis <4 || dis == 7){
    forwardjet = LorentzVector(0,0,0,0); 
    string unusedjets = result.get_unusedJets();
    //cout<<"Not used jets "<<unusedjets<<endl;
    for(unsigned int i= 0; i<event.jets->size(); i++){
      //cout<<!bool(unusedjets[i])<<" eta jet "<< event.jets->at(i).eta()<<endl;
      if(7>i){
	if(!unusedjets[i]){
	   forwardnumber++;
	   if(fabs(event.jets->at(i).eta())>=2)forwardnumber_eta2++;
        }
	if(!unusedjets[i] && fabs(event.jets->at(i).eta())>fabs(forwardjet.eta()))
	  forwardjet = event.jets->at(i).v4();
      }
      else if(fabs(event.jets->at(i).eta())>fabs(forwardjet.eta())){
	forwardjet = event.jets->at(i).v4();
      	forwardnumber++;
	if(fabs(event.jets->at(i).eta())>=2)forwardnumber_eta2++;
      }
      else{
	forwardnumber++;
	if(fabs(event.jets->at(i).eta())>=2)forwardnumber_eta2++;
      }
    }
    
  }
  else if(dis == 6){
    forwardjet = LorentzVector(0,0,0,0); 
    LorentzVector wtagjet = result.get_wHad();
    LorentzVector topjet = result.get_topJets();
    for(unsigned int i= 0; i<event.jets->size(); i++){
      if(deltaR(wtagjet,event.jets->at(i).v4())>0.8 && deltaR(topjet,event.jets->at(i).v4())> 0){
	forwardnumber++;
	if(fabs(event.jets->at(i).eta())>=2)forwardnumber_eta2++;
	if(fabs(event.jets->at(i).eta())>fabs(forwardjet.eta()) )
	  forwardjet = event.jets->at(i).v4();
      }
    }
  }
  else if(dis == 5 || dis == 4 || dis == 8){
    forwardjet = LorentzVector(0,0,0,0); 
    LorentzVector toptagjet = result.get_topHad();
     for(unsigned int i= 0; i<event.jets->size(); i++){
       if(deltaR(toptagjet,event.jets->at(i).v4())>0.8 || fabs(event.jets->at(i).eta())>2.4){
	 forwardnumber++;
	 if(fabs(event.jets->at(i).eta())>=2)forwardnumber_eta2++;
	 if(fabs(event.jets->at(i).eta())>fabs(forwardjet.eta())){
	   forwardjet = event.jets->at(i).v4();
	 }
       }
     }
  }
  result.set_num_forwardjets(forwardnumber);
  result.set_num_forwardjets_eta2(forwardnumber_eta2);
  result.set_forwardJet(forwardjet);
  double jetiso = -1;
  if(forwardjet.E()>0)
    for(auto jet : *event.jets)
      if((jetiso > deltaR(jet.v4(),forwardjet) || jetiso ==-1) && deltaR(jet.v4(),forwardjet)>0)jetiso = deltaR(jet.v4(),forwardjet);
  result.set_jetIso(jetiso);

  /*/
  if(dis == 1){
      if(result.get_Mass()==-1){
	cout<<"No reconstructed mass, event has number of jets: "<<event.jets->size()<<endl;	
	for(auto j : *event.jets)
	  cout<<j.pt()<<" "<<j.eta()<<" "<<j.phi()<<endl;
	assert(1==0);
      }
  }
  /*/
  
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
  
  double mean_topLep = 170.1, sigma_topLep = 19.1;
  double mean_topHad = 172.6, sigma_topHad = 14.29;
  boost::math::chi_squared two_param(1);

  for(auto hyp :  event.get(hyps)){
    const LorentzVector & whad    = hyp.get_wHad();
    const LorentzVector & wlep    = hyp.get_wLep();
    const LorentzVector & topJets = hyp.get_topJets();
    
    //this is the chi2
    double ttbar = pow(((topJets+wlep).M()-mean_topLep)/sigma_topLep,2)+pow((whad.M()-mean_topHad)/sigma_topHad,2);
    ttbar =  boost::math::cdf(two_param,ttbar);
    if(ttbar > ttbarchi || ttbarchi == -1){
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
  LorentzVector bprime(0,0,0,0);
  if(dis==7) genPart = event.get(gen);

  //cout<<"================================================"<<endl;
  // values for phi+eta and after // for deltaR
  boost::math::chi_squared four_param(4); //3
  boost::math::chi_squared three_param(3); //2
  boost::math::chi_squared two_param(2); //1
 
  for(auto hyp :  event.get(hyps)){
    const LorentzVector & whad    = hyp.get_wHad();
    const LorentzVector & wlep    = hyp.get_wLep();
    const LorentzVector & topJets = hyp.get_topJets();

    //cout<<"***************************"<<endl;
    //cout<<" whad "<<whad.pt()<<" wlep "<<wlep.pt()<<" topjet "<<topJets.pt()<<endl;
    //cout<< "bprime chi "<<bprimechi<<endl;

      
    double lepTop =999999999;
    double hadTop =999999999;

    //double mean_topLep = 175., sigma_topLep = 18.7, mean_wHad = 81., sigma_wHad = 9.6; 
    //double mean_topHad = 167., sigma_topHad = 21.93;
    double mean_topLep = 170.1,  sigma_topLep = 19.1, mean_wHad = 85.5, sigma_wHad = 8.7; 
    double mean_topHad = 172.6,  sigma_topHad = 14.29;
    double mean_distance = 3.14, sigma_distance = 0.8;
    double mean_ptratio = 1, sigma_ptratio = 0.5;
    double mean_phi = 3.14,  sigma_phi = 0.8;
    double mean_eta = 0.0,   sigma_eta = 0.8;
    
    double factor = 1.;

    //double distance_tlep = pow(deltaR(whad,wlep+topJets)-mean_distance,2)/(sigma_distance*sigma_distance);
    //double distance_thad = pow(deltaR(wlep,whad+topJets)-mean_distance,2)/(sigma_distance*sigma_distance);
    double distance_tlep = pow(deltaPhi(whad,wlep+topJets)-mean_phi,2)/(sigma_phi*sigma_phi)+pow(fabs(whad.eta()-(wlep+topJets).eta())-mean_eta,2)/(sigma_eta*sigma_eta);
    double distance_thad = pow(deltaPhi(wlep,whad+topJets)-mean_phi,2)/(sigma_phi*sigma_phi)+pow(fabs(wlep.eta()-(whad+topJets).eta())-mean_eta,2)/(sigma_eta*sigma_eta);

    
    sigma_topLep *= factor; sigma_wHad *= factor; sigma_topHad *=factor; sigma_distance *= factor; sigma_ptratio *=factor;

    if(topJets.pt()>0&&( dis==1 || dis ==2))lepTop = (((topJets+wlep).M()-mean_topLep)*((topJets+wlep).M()-mean_topLep)/(sigma_topLep*sigma_topLep)+(whad.M()-mean_wHad)*(whad.M()-mean_wHad)/(sigma_wHad*sigma_wHad)+distance_tlep+pow(whad.pt()/(wlep+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio));//*0.25*0.25;

    if(dis==1 || dis ==3)hadTop = (((topJets+whad).M()-mean_topHad)*((topJets+whad).M()-mean_topHad)/(sigma_topHad*sigma_topHad)+distance_thad+pow(wlep.pt()/(whad+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio));///9;
    
    lepTop = boost::math::cdf(four_param,lepTop);
    hadTop = boost::math::cdf(three_param,hadTop);
    
    if(hyp.get_whad_num()==1){
      if(topJets.pt()>0.){	
        if(dis==1 || dis ==2){
	  lepTop = (((topJets+wlep).M()-mean_topLep)*((topJets+wlep).M()-mean_topLep)/(sigma_topLep*sigma_topLep)+ distance_tlep+pow(whad.pt()/(wlep+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio));//9;
	  lepTop = boost::math::cdf(three_param,lepTop);
        }
        if(dis==1 || dis ==3){
	  // Andrews request
	  //if(hyp.get_top_num()==1 && (topJets+whad).pt()<400) continue;
	  hadTop = (((topJets+whad).M()-mean_topHad)*((topJets+whad).M()-mean_topHad)/(sigma_topHad*sigma_topHad)+distance_thad+pow(wlep.pt()/(whad+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio));//9; 
  	  hadTop = boost::math::cdf(three_param,hadTop);
        }
     }
      else{
	// Andrews request 
	//if(whad.pt()<800)continue; 
	if(dis==1 || dis ==3){
	  hadTop = (distance_thad+pow(wlep.pt()/(whad+topJets).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio));//4; 
	  hadTop = boost::math::cdf(two_param,hadTop);
	}
      }
   }     
   
     //cout<<"chi2 p top lep "<<lepTop<<endl;//" p-value "<< boost::math::cdf(four_param,lepTop)<<endl;
     //cout<<"chi2 p top had "<<lepTop<<endl;//" p-value "<< boost::math::cdf(three_param,lepTop)<<endl;

    /*/
    if(deltaR(whad+topJets,wlep)<2.4)
       hadTop +=10000000;
    if(deltaR(wlep+topJets,whad)<2.4)
       lepTop +=10000000;
    /*/
    //cout<<" chi top lep "<<lepTop<<" top had "<<hadTop<<endl;
    
    if(dis==1){
      // Top lep == 11 / Top had == 12
      double combochi = lepTop > hadTop ? hadTop : lepTop;
      int recoType_help = lepTop > hadTop ? 12 : 11;
      if(combochi < bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = combochi;
	//cout<<" combo chi "<<combochi<<" stored chi "<< bprimechi<<endl;
	recoType = recoType_help;
	mass = sqrt((wlep+whad+topJets).M2());
	bprime = wlep+whad+topJets;
      }
    }
    else if(dis==2){
      if(lepTop<bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = lepTop;
	recoType = 11;
	mass = sqrt((wlep+whad+topJets).M2());
	bprime =wlep+whad+topJets;
      }
    }
    else if(dis==3){
      if(hadTop<bprimechi || bprimechi ==-1){
	bestHyp = hyp;
	bprimechi = hadTop;
	recoType = 12;
	mass = sqrt((wlep+whad+topJets).M2());
	bprime = wlep+whad+topJets;
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
	gen_chi = pow(deltaR(genPart.get_topHad(),(whad+topJets)),2)/0.25+pow((whad+topJets).M()-genPart.get_topHad().M(),2)/1;
	if(deltaR(genPart.get_topHad(),(whad+topJets))>0.4  || deltaR(genPart.get_wHad(),whad)>0.4 || deltaR(genPart.get_wLep(),wlep)>0.4|| round((genPart.get_topHad()+genPart.get_wLep()).M()) != round(genPart.get_bprime().M()))
	  gen_chi +=99999;
	//cout<<"Gen B Mass from TW "<<round((genPart.get_topHad()+genPart.get_wLep()).M())<<" from Gen particle "<<round(genPart.get_bprime().M())<<" "<< bool(round((genPart.get_topHad()+genPart.get_wLep()).M()) != round(genPart.get_bprime().M())) <<endl;      
      }
      else{ 
	gen_chi = pow(deltaR(genPart.get_topLep(),wlep+topJets),2)/0.16+pow(deltaR(genPart.get_wHad(),whad),2)/0.16+pow((wlep+topJets).M()-genPart.get_topLep().M(),2)/8+pow(whad.M()-genPart.get_wHad().M(),2)/8;
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
  bestHyp.set_bprime(bprime);
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
  double subjetbtag = -1;
  LorentzVector bprime(0,0,0,0);

  double mean_topHad = 172.6, sigma_topHad = 14.29;
  double mean_distance = 3.14, sigma_distance = 0.8;
  double mean_ptratio = 1, sigma_ptratio = 0.4;
  boost::math::chi_squared three_param(2);
  
  for(auto hyp :  event.get(hyps)){
    const LorentzVector & wlep = hyp.get_wLep();
    const LorentzVector & topHad = hyp.get_topHad();


    // chi2/nodf
    double chi2 = pow((topHad.M()-mean_topHad)/sigma_topHad,2)+pow((deltaR(wlep,topHad)-mean_distance)/sigma_distance,2)+pow((wlep.pt()/topHad.pt()-mean_ptratio)/sigma_ptratio,2);
    double prob = boost::math::cdf(three_param,chi2);
    
    if(prob>bprimechi || bprimechi==-1){
      bestHyp=hyp;
      subjetbtag = hyp.get_btag_discriminator();
      bprimechi=prob;
      recoType=2;
      mass = sqrt((wlep+topHad).M2());
      bprime = wlep+topHad;
    }
  }
  if(bprimechi==-1) return bestHyp;
  bestHyp.set_btag_discriminator(subjetbtag);
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_Mass(mass);
  bestHyp.set_bprime(bprime);
  return bestHyp;
}

BprimeContainer BprimeDiscriminator::wTag_dis(uhh2::Event & event){
  BprimeContainer bestHyp;
  double bprimechi =-1;
  int recoType = -1;
  double mass =-1;
  LorentzVector bprime(0,0,0,0);
  LorentzVector empty(0,0,0,0);
 
  
  for(auto hyp :  event.get(hyps)){
    const LorentzVector & wlep = hyp.get_wLep();
    const LorentzVector & wHad = hyp.get_wHad();
    if(wHad.pt()==0)continue;
    
    double chilepTop=0;
    double chihadTop=0;
    
    //double mean_topHad = 167., sigma_topHad = 21.93;
    double mean_topLep = 170.1, sigma_topLep = 19.1, mean_wHad = 85.5, sigma_wHad = 8.7; 
    double mean_topHad = 172.6, sigma_topHad = 14.29;
    double mean_distance = 3.14, sigma_distance = 0.15;
    double mean_ptratio = 1, sigma_ptratio = 0.5;

    for(unsigned int i =0; i<event.jets->size()+1;i++){
      LorentzVector toplep = wlep;
      LorentzVector tophad = wHad;	
      LorentzVector jet(0,0,0,0);
      if(i!=event.jets->size()){
	jet = (*event.jets)[i].v4();
	if(abs(jet.eta())>2.4 || deltaR(jet,wHad)<=1.2 || jet.pt()<30)
	//if(abs(jet.eta())>2.4)
	  continue;
	toplep += jet;
	tophad += jet;
      }
 
      chilepTop = ((toplep.M()-mean_topLep)*(toplep.M()-mean_topLep)/(sigma_topLep*sigma_topLep)+(wHad.M()-mean_wHad)*(wHad.M()-mean_wHad)/(sigma_wHad*sigma_wHad)+pow(deltaR(wHad,toplep)-mean_distance,2)/(sigma_distance*sigma_distance)+pow(wHad.pt()/(toplep).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio))*0.25*0.25;
      chihadTop = ((tophad.M()-mean_topHad)*(tophad.M()-mean_topHad)/(sigma_topHad*sigma_topHad)+pow(deltaR(wlep,tophad)-mean_distance,2)/(sigma_distance*sigma_distance)+pow(wlep.pt()/(tophad).pt()-mean_ptratio,2)/(sigma_ptratio*sigma_ptratio))/9; 


      if(chilepTop < chihadTop && i==event.jets->size())
	continue;

      if(chilepTop<bprimechi || bprimechi==-1){
	bestHyp=hyp;
	bestHyp.set_topJets(jet);
	bestHyp.set_topLep(toplep);
	bprimechi=chilepTop;
	  recoType=61;
	  mass = sqrt((wHad+toplep).M2());
	  bprime = wHad+toplep;
      }
      if(chihadTop<bprimechi || bprimechi==-1){
	bestHyp=hyp;
	bestHyp.set_topJets(jet);
	bestHyp.set_topHad(tophad);
	bprimechi=chihadTop;
	recoType=62;
	mass = sqrt((wlep+tophad).M2());
	bprime = wlep+tophad;
      }
    }
  }
  if(bprimechi==-1) return bestHyp;
  bestHyp.set_chiVal(bprimechi);
  bestHyp.set_RecoTyp(recoType);
  bestHyp.set_Mass(mass);
  bestHyp.set_bprime(bprime);
  return bestHyp;
}
