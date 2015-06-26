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
  vector<BprimeContainer> recoWHyps;
  vector<BprimeContainer> recoHyps;
  BprimeContainer tmp_hyp;
  //reconstruct leptonic W
  LorentzVector lep(0,0,0,0);
  lep = event.get(primlep).v4();  
  //cout<<"lep pt "<<lep.pt()<<endl;
  vector<LorentzVector> neutrinos = NeutrinoReconstruction(lep, event.met->v4());
  //reconstruct hadronic W and add leptonic W
  for(unsigned int K=1; K<jets.size(); ++K){
    unsigned int N = jets.size();
    string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's 
    //permute bitmask
    do {
      LorentzVector whad(0,0,0,0);
      for (unsigned int i = 0; i < N; ++i){ // [0..N-1] integers
	if(bitmask[i]) whad = whad+jets.at(i).v4();
	if(whad.M()>600) break;
      }
      for(auto & neutrino:neutrinos){
	if(whad.M()<400 && whad.M()>40){
	  tmp_hyp.set_wJets(bitmask);
	  tmp_hyp.set_wHad(whad);
	  tmp_hyp.set_wLep(neutrino+lep);
	  recoWHyps.emplace_back(tmp_hyp);
	}
      }
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  }
  //cout<<"Size of reco hyps "<<recoWHyps.size()<<endl;

  //permutation of top jets
  for(auto & hyp : recoWHyps){
    recoHyps.emplace_back(hyp);//hypothersis without any jets!
    vector<Jet> unusedJets;
    //cout<<hyp.get_wJets().size()<<endl;
    for(unsigned int i=0; i<jets.size();++i){
      if(!hyp.get_wJets()[i]) unusedJets.emplace_back(jets.at(i));
    }
    for(unsigned int K=1; K<unusedJets.size(); ++K){
      unsigned int N = unusedJets.size();
      string bitmask(K, 1); // K leading 1's
      bitmask.resize(N, 0); // N-K trailing 0's 
      //permute bitmask
      do {
	LorentzVector top(0,0,0,0);
	for (unsigned int i = 0; i < N; ++i){ // [0..N-1] integers
	  if (bitmask[i]) top = top+unusedJets.at(i).v4();
	  //if (bitmask[i]) cout <<i<<endl;
	}
	if((top+hyp.get_wHad()).M()>400 && (top+hyp.get_wLep()).M()>400) continue;
	hyp.set_topJets(top);
	//hyp.set_topJets(bitmask); //intendet for a string to see which jets where associated
  	recoHyps.emplace_back(hyp);
      } while(std::prev_permutation(bitmask.begin(), bitmask.end()));   
    }
  }
  if(recoHyps.size()==0) return false;
  event.set(hypothesis, move(recoHyps));
  return true;
}

void BprimeReco::comb(int N, int K)
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
  vector<LorentzVector> neutrinos = NeutrinoReconstruction(lep, event.met->v4());
  //vector<Jet>* jets = event.jets;
  vector<BprimeContainer> recoHyps;
  BprimeContainer tmp_hyp;

  if(topjets.size()<1 || neutrinos.size()<1)return false;
  for(auto & topjet : topjets){
    tmp_hyp.set_topHad(topjet.v4());
    for(auto & neutrino : neutrinos){
      if(deltaR(topjet.v4(),neutrino+lep)>dRmin){
	tmp_hyp.set_wLep(neutrino+lep);
	recoHyps.push_back(tmp_hyp);
      }
    }
  }
  if(recoHyps.size()<1)return false;
  event.set(hypothesis, move(recoHyps));
  return true;
}


template<typename T>
  bool BprimeReco::passes_id(const T & object, const Event & event, const boost::optional<std::function<bool (const T &, const Event & )>> & object_id){
  if((*object_id)(object, event))
    return true;
  return false;
}
