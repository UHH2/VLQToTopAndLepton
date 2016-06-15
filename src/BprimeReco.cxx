#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"

#include "UHH2/common/include/TTbarReconstruction.h"


#include <algorithm>
#include <iostream>
#include <string>

using namespace std;
using namespace uhh2;

BprimeReco::BprimeReco(uhh2::Context & ctx, const std::string & label){
  hypothesis = ctx.declare_event_output<vector<BprimeContainer>>(label);
  primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  topjetCollBool = false;
  jetCollBool =false;
  jetId=boost::none;
  topjetId=boost::none;
  FitNeutrino = NeutrinoFit();
}

bool BprimeReco::process(uhh2::Event & event){
  cout<<"Don't use the 'process' function for the event reconstruction for Bprime"<<endl;
  return false;
}

bool BprimeReco::set_topjetCollection(uhh2::Context & ctx,const string & topjetCollectionName){
  topjet_collection = ctx.get_handle<vector<TopJet>>(topjetCollectionName);
  topjetCollBool =true;
  return true;
}
  
bool BprimeReco::set_jetCollection(uhh2::Context & ctx,const string & jetCollectionName){
 jet_collection = ctx.get_handle<vector<Jet>>(jetCollectionName);
 jetCollBool =true;
 return true;
}

bool BprimeReco::massReco(uhh2::Event & event){
  vector<Jet> jets = jetCollBool ? event.get(jet_collection) : *event.jets;
  if(jetId){
    vector<Jet> selected_jets; 
    for(unsigned int i =0; i<jets.size(); i++){
	if(passes_id(jets.at(i),event,jetId)) selected_jets.push_back(jets.at(i));
    }
    selected_jets.swap(jets);
  }
  //vector<BprimeContainer> recoWHyps;
  vector<BprimeContainer> recoHyps;
  LorentzVector lep = event.get(primlep).v4();  
  //vector<LorentzVector> wleps = NeutrinoReconstruction(lep, event.met->v4());
  vector<LorentzVector> wleps = FitNeutrino.NeutrinoFitPolar(lep,event.met->v4());
  for(unsigned int i=0; i<wleps.size();++i) wleps[i] +=  lep;
  vector<BprimeContainer> recoWHyps = reconstruct_WHyps(jets,wleps);
  //permutation of top jets
  for(auto & hyp : recoWHyps){
    string jet_string = hyp.get_wJets();
    hyp.set_unusedJets(jet_string);//setting unused Jets to wJets for hyp with no Jets!
    recoHyps.emplace_back(hyp);//hypothersis without any additional jets!
    vector<Jet> unusedJets;
    vector<double> unusedjetsMap;
    for(unsigned int i=0; i<jets.size();++i){
      if(!jet_string[i]){
	unusedJets.emplace_back(jets.at(i));
	unusedjetsMap.emplace_back(i);
      }
    }  
    for(unsigned int K=1; K<unusedJets.size(); ++K){
      string hyp_jet_string =  jet_string;
      unsigned int N = unusedJets.size();
      string bitmask(K, 1); // K leading 1's
      bitmask.resize(N, 0); // N-K trailing 0's 
      vector<LorentzVector> jets_top;
      //permute bitmask
      do {
	LorentzVector top(0,0,0,0);
	for (unsigned int i = 0; i < N; ++i){ // [0..N-1] integers
	  if (bitmask[i]){
	    top = top+unusedJets.at(i).v4();
	    jets_top.push_back(unusedJets.at(i).v4());
	  }
	  //if (bitmask[i]) cout <<i<<endl;
	}
	if((top+hyp.get_wHad()).M()>500 && (top+hyp.get_wLep()).M()>500) continue;
	hyp.set_topJets(top);
	//hyp.set_topLorentz(jets_top);
	for(unsigned int it = 0; it<bitmask.size(); it++){
	  if(bitmask[it]){
	    hyp_jet_string[unusedjetsMap[it]]=1;
	  }
	}
	//cout<<"wJets allJets size "<< jets.size()<<endl;
	//for(unsigned int ip =0; ip< hyp.get_wJets().size(); ip++)
	//  cout<<bool(hyp.get_wJets()[ip])<<" "<<bool(hyp_jet_string[ip])<<endl;
	//cout<<"==================="<<endl;
	hyp.set_unusedJets(hyp_jet_string);
  	recoHyps.emplace_back(hyp);
      } while(std::prev_permutation(bitmask.begin(), bitmask.end()));   
    }
  }
  int size_recoHyps = recoHyps.size();
  event.set(hypothesis, move(recoHyps));
  return size_recoHyps > 0 ? true : false;
}
//just an example not used anywhere
void BprimeReco::comb(int N, int K)//prototyp of the permutation algorithem
{
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
 
    // print integers and permute bitmask
    do {
      for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
	  if (bitmask[i]) cout<<" "<<i;
        }
      std::cout << std::endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}


bool BprimeReco::TopJetReco(Event & event, double dRmin){
  vector<TopJet> topjets = topjetCollBool ?  event.get(topjet_collection) : *event.topjets;
  if(topjetId){
    vector<TopJet> selected_topjets;
    for(unsigned int i =0; i<topjets.size(); i++){
      if(passes_id(topjets.at(i),event,topjetId)) selected_topjets.push_back(topjets.at(i));
    }
    selected_topjets.swap(topjets);
  }	  
  LorentzVector lep = event.get(primlep).v4();  
  //vector<LorentzVector> neutrinos = NeutrinoReconstruction(lep, event.met->v4());
  vector<LorentzVector> neutrinos = FitNeutrino.NeutrinoFitPolar(lep,event.met->v4());

  vector<BprimeContainer> recoHyps;
  BprimeContainer tmp_hyp;
  //if(topjets.size()==0||neutrinos.size()==0) return false;
  for(auto & topjet : topjets){
    const vector<Jet> & Top_subjets = topjet.subjets();
    LorentzVector TopCand(0,0,0,0);
    for(auto &subjet :Top_subjets)
      TopCand += subjet.v4();
    tmp_hyp.set_topHad(TopCand);
    //tmp_hyp.set_topJets(TopCand);
    for(auto & neutrino : neutrinos){
      if(deltaR(topjet.v4(),neutrino+lep)>dRmin){
	tmp_hyp.set_wLep(neutrino+lep);
	recoHyps.push_back(tmp_hyp);
      }
    }
  }
  int size_recoHyps = recoHyps.size();
  event.set(hypothesis,move(recoHyps));
  return size_recoHyps > 0 ? true : false;
}

bool BprimeReco::BTagReco(Event & event){
  if(!jetId){
    cerr<<"You have to set a jet Id which to veto"<<endl;
    assert(0==1);
  }
  vector<Jet>* jets = jetCollBool ? &event.get(jet_collection) : event.jets;
  vector<Jet> selected_jets; 
  vector<Jet> bjets; 
  for(unsigned int i =0; i<jets->size(); i++){
    if(!passes_id(jets->at(i),event,jetId)) selected_jets.push_back(jets->at(i));
    else bjets.push_back(jets->at(i));
  }
  LorentzVector lep = event.get(primlep).v4();  
  //vector<LorentzVector> wleps = NeutrinoReconstruction(lep, event.met->v4());
  vector<LorentzVector> wleps = FitNeutrino.NeutrinoFitPolar(lep,event.met->v4());
  for(unsigned int i=0; i<wleps.size();++i) wleps[i] +=  lep;
  vector<BprimeContainer> recoWHyps = reconstruct_WHyps(selected_jets,wleps);
  vector<BprimeContainer> recoHyps;
  for(auto & hyp : recoWHyps){
    for(auto & bjet :  bjets){
      hyp.set_topJets(bjet.v4());
      recoHyps.push_back(hyp);
    }
  }
  int size_recoHyps = recoHyps.size();
  event.set(hypothesis, move(recoHyps));
  return size_recoHyps > 0 ? true : false;
}

template<typename T>
  bool BprimeReco::passes_id(const T & object, const Event & event, const boost::optional<std::function<bool (const T &, const Event & )>> & object_id){
  if((*object_id)(object, event))
    return true;
  return false;
}

vector<BprimeContainer> BprimeReco::reconstruct_WHyps(const std::vector<Jet> & jets, const std::vector<LorentzVector> & Wleps, double cutoff_WHad_min, double cutoff_WHad_max){
  vector<BprimeContainer> recoWHyps;
  BprimeContainer tmp_hyp;
  //reconstruct hadronic W and add leptonic W
  for(unsigned int K=1; K<jets.size(); ++K){
    unsigned int N = jets.size();
    string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's 
    vector<LorentzVector> whad_jets;
    //permute bitmask
    do {
      LorentzVector whad(0,0,0,0);
      for (unsigned int i = 0; i < N; ++i){ // [0..N-1] integers
	if(bitmask[i]){ 
	  whad = whad+jets.at(i).v4();
	  whad_jets.push_back(jets.at(i).v4());
	}
	if(whad.M()>cutoff_WHad_max*1.25) break;
      }
      for(auto & wlep:Wleps){
	if(whad.M()<cutoff_WHad_max && whad.M()>cutoff_WHad_min){
	  tmp_hyp.set_wJets(bitmask);
	  tmp_hyp.set_wHad(whad);
	  //tmp_hyp.set_wHadJets(whad_jets);
	  tmp_hyp.set_wLep(wlep);
	  recoWHyps.emplace_back(tmp_hyp);
	}
      }
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  }
  return recoWHyps;
}

bool BprimeReco::hadronicW(uhh2::Event & event, double dRmin){
  vector<TopJet> wjets = topjetCollBool ?  event.get(topjet_collection) : *event.topjets;
  vector<Jet> jets = jetCollBool ? event.get(jet_collection) : *event.jets;
  if(wjetId){
    vector<TopJet> selected_wjets;
    for(unsigned int i =0; i<wjets.size(); i++){
      if(passes_id(wjets.at(i),event,wjetId)) selected_wjets.push_back(wjets.at(i));
    }
    selected_wjets.swap(wjets);
  }	  
  LorentzVector lep = event.get(primlep).v4();  
  //vector<LorentzVector> neutrinos = NeutrinoReconstruction(lep, event.met->v4());
  vector<LorentzVector> neutrinos = FitNeutrino.NeutrinoFitPolar(lep,event.met->v4());
  vector<BprimeContainer> recoHyps;
  BprimeContainer tmp_hyp;
  for(auto & wjet : wjets){
    const vector<Jet> & w_subjets = wjet.subjets();
    LorentzVector WCand(0,0,0,0);
    for(auto &subjet : w_subjets)
      WCand += subjet.v4();
    tmp_hyp.set_wHad(WCand);
    //tmp_hyp.set_topJets(TopCand);
    for(auto & neutrino : neutrinos){
      if(deltaR(wjet.v4(),neutrino+lep)>dRmin){
	tmp_hyp.set_wLep(neutrino+lep);
	for(auto & jet :jets){
	  if(deltaR(jet.v4(),WCand)>dRmin){
	    tmp_hyp.set_topJets(jet.v4());
	    tmp_hyp.set_topLep(jet.v4()+neutrino+lep);
	    recoHyps.push_back(tmp_hyp);
	  }
	}
      }
    }
  }
  int size_recoHyps = recoHyps.size();
  event.set(hypothesis,move(recoHyps));
  return size_recoHyps > 0 ? true : false;
}
