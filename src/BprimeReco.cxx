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
}

bool BprimeReco::process(uhh2::Event & event){
  cout<<"Don't use the 'process' function for the event reconstruction for Bprime"<<endl;
  return false;
}

bool BprimeReco::massReco(uhh2::Event & event){
  vector<Jet>* jets = event.jets;
  vector<BprimeContainer> recoWHyps;
  vector<BprimeContainer> recoHyps;
  BprimeContainer tmp_hyp;

  //reconstruct leptonic W
  LorentzVector lep(0,0,0,0);
  lep = event.get(primlep).v4();  
  //cout<<"lep pt "<<lep.pt()<<endl;
  vector<LorentzVector> neutrinos = NeutrinoReconstruction(lep, event.met->v4());

  //reconstruct hadronic W and add leptonic W
  for(unsigned int K=1; K<jets->size(); ++K){
    unsigned int N = jets->size();
    string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's 
    //permute bitmask
    do {
      LorentzVector whad(0,0,0,0);
      for (unsigned int i = 0; i < N; ++i){ // [0..N-1] integers
	  if (bitmask[i]) whad = whad+jets->at(i).v4();
	  if(whad.M()>600) break;
      }
      for(auto & neutrino:neutrinos){
	if(whad.M()<500 && whad.M()>40 ){
	  tmp_hyp.set_wJets(bitmask);
	  tmp_hyp.set_wHad(whad);
	  tmp_hyp.set_wLep(neutrino+lep);
	  recoWHyps.emplace_back(tmp_hyp);
	  //cout<<whad.M()<<" "<<(neutrino+lep).M()<<endl;
	}
      }
      //cout<<whad.M()<<endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  }

  //permutation of top jets
  for(auto & hyp : recoWHyps){
    vector<Jet> unusedJets;
    for(unsigned int i=0; i<jets->size();++i)
      if(!hyp.get_wJets()[i]) unusedJets.emplace_back(jets->at(i));
    for(unsigned int K=1; K<unusedJets.size(); ++K){
      unsigned int N = unusedJets.size();
      string bitmask(K, 1); // K leading 1's
      bitmask.resize(N, 0); // N-K trailing 0's 
      //permute bitmask
      do {
	LorentzVector top(0,0,0,0);
	for (unsigned int i = 0; i < N; ++i) // [0..N-1] integers
	  if (bitmask[i]) top = top+unusedJets.at(i).v4();
	if((top+hyp.get_wHad()).M()>400 && (top+hyp.get_wLep()).M()>400) break;
	hyp.set_topJets(top);
	//hyp.set_topJets(bitmask); //intendet for a string to see which where taken
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
