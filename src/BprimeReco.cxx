#include "UHH2/VLQToTopAndLepton/include/BprimeReco.h"

#include "UHH2/common/include/TTbarReconstruction.h"
#include <boost/algorithm/string.hpp>    

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
  bjetId = CSVBTag(CSVBTag::WP_MEDIUM);

  string helper_lable = label;
  boost::algorithm::to_lower(helper_lable);
  if(boost::algorithm::contains(helper_lable,"wtag")){
    wtag =true;
    file = TFile::Open( "/nfs/dust/cms/user/gonvaq/CMSSW/CMSSW_8_0_24_patch1/src/UHH2/VLQToTopAndLepton/data/puppiCorr.root","READ");
    puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
    puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
    puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");
    file->Close();
  }
}

BprimeReco::~BprimeReco(){
  if(wtag)
    file->Close();
}

int BprimeReco::count_bjets(vector<Jet> & jets, Event & event){
  int counter = 0;
  for(unsigned int i =0; i<jets.size(); i++){
    if(passes_id(jets.at(i),event,bjetId))
      counter++;
  }
  return counter;
}

double BprimeReco::forwardeta(vector<Jet> & jets, Event & event){
  double eta = 0;
  if(forwardJetId){
    for(unsigned int i =0; i<jets.size(); i++){
      if(passes_id(jets.at(i),event,forwardJetId))
	if(eta<fabs(jets.at(i).eta()))eta =fabs(jets.at(i).eta());
    }
  }
  else{
    for(unsigned int i =0; i<jets.size(); i++){
      if(eta<fabs(jets.at(i).eta()))eta =fabs(jets.at(i).eta());
    }
  }
  return eta;
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
  int numberOfBjets = count_bjets(*event.jets,event);
  double forwardjeteta = forwardeta(*event.jets,event); 

  vector<Jet> jets = jetCollBool ? event.get(jet_collection) : *event.jets;
  if(jetId){
    vector<Jet> selected_jets; 
    for(unsigned int i =0; i<jets.size(); i++){
	if(passes_id(jets.at(i),event,jetId)) selected_jets.push_back(jets.at(i));
    }
    selected_jets.swap(jets);
  }
  //if(jets.size()>12)jets.resize(12);
  if(jets.size()>7)jets.resize(7);

  
  //vector<BprimeContainer> recoWHyps;
  vector<BprimeContainer> recoHyps;
  LorentzVector lep = event.get(primlep).v4();  
  //vector<LorentzVector> wleps = NeutrinoReconstruction(lep, event.met->v4());
  vector<LorentzVector> wleps = FitNeutrino.NeutrinoFitPolar(lep,event.met->v4());
  for(unsigned int i=0; i<wleps.size();++i) wleps[i] +=  lep;
  vector<BprimeContainer> recoWHyps = reconstruct_WHyps(jets,wleps);
  //permutation of top jets
  for(auto & hyp : recoWHyps){
    hyp.set_forwardJetEta(forwardjeteta);
    hyp.set_EventBtagNumber(numberOfBjets);
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
      double bdiscrimi = -1;
      vector<LorentzVector> jets_top;
      //permute bitmask
      do {
        hyp.set_num_top(0);
	LorentzVector top(0,0,0,0);
	for (unsigned int i = 0; i < N; ++i){ // [0..N-1] integers
	  if (bitmask[i]){
	    top = top+unusedJets.at(i).v4();
            hyp.add_num_top();
	    jets_top.push_back(unusedJets.at(i).v4());
	    bdiscrimi = unusedJets.at(i).btag_combinedSecondaryVertex();
	  }
	  //if((top+hyp.get_wHad()).M()>500 && (hyp.get_top_num()+hyp.get_whad_num())>1)
	  //  break;
	  //if (bitmask[i]) cout <<i<<endl;
	}
	//if((top+hyp.get_wHad()).M()>500 && (top+hyp.get_wLep()).M()<50) continue;
	hyp.set_topJets(top);
	if(hyp.get_btag_discriminator() < bdiscrimi)hyp.set_btag_discriminator(bdiscrimi);
	for(unsigned int it = 0; it<bitmask.size(); it++){
	  if(bitmask[it]){
	    hyp_jet_string[unusedjetsMap[it]]=1;
	  }
	}
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
  int numberOfBjets = count_bjets(*event.jets,event);
  double forwardjeteta = forwardeta(*event.jets,event); 

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
  tmp_hyp.add_num_top();
  tmp_hyp.set_forwardJetEta(forwardjeteta);
  tmp_hyp.set_EventBtagNumber(numberOfBjets);
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
  //cout<<"reco hypothesis size "<<size_recoHyps<<endl;
  event.set(hypothesis,move(recoHyps));
  return size_recoHyps > 0 ? true : false;
}

bool BprimeReco::BTagReco(Event & event){
  int numberOfBjets = count_bjets(*event.jets,event);
  double forwardjeteta = forwardeta(*event.jets,event); 

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
    hyp.set_forwardJetEta(forwardjeteta);
    hyp.set_EventBtagNumber(numberOfBjets);
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
  //cout<<"======================="<<endl;
  for(unsigned int K=1; K<jets.size()+1; ++K){
    double bdiscrimi = -1;
    unsigned int N = jets.size();
    string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's 
    vector<LorentzVector> whad_jets;
    //permute bitmask
    do{
      tmp_hyp.set_num_whad(0);
      LorentzVector whad(0,0,0,0);
      for (unsigned int i = 0; i < N; ++i){ // [0..N-1] integers
	if(bitmask[i]){
	  //cout<<1<<" ";
	  whad = whad+jets.at(i).v4();
          tmp_hyp.add_num_whad();
	  whad_jets.push_back(jets.at(i).v4());
	  bdiscrimi = jets.at(i).btag_combinedSecondaryVertex();	   
	}
	//else
	//  cout<<0<<" ";
	//if(whad.M()>cutoff_WHad_max*1.25 || i > 4) break;
      }
      //cout<<endl;
      //if(tmp_hyp.get_whad_num()<2)continue;
      for(auto & wlep : Wleps){
	if(whad.M()<cutoff_WHad_max && whad.M()>cutoff_WHad_min){
	  tmp_hyp.set_wJets(bitmask);
	  tmp_hyp.set_wHad(whad);
	  if(tmp_hyp.get_btag_whad() < bdiscrimi)tmp_hyp.set_btag_whad(bdiscrimi);
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
    if(deltaR(wjet.v4(),lep)<= dRmin) continue;
    LorentzVector raw_v4(0,0,0,0);
    for(auto sub :wjet.subjets())
      raw_v4 += sub.v4()*sub.JEC_factor_raw();

    double mass = raw_v4.M()*calc_wtag_corr(wjet.eta(), wjet.pt());
    wjet.set_softdropmass(mass);
    //cout<<"softdrop mass raw "<<raw_v4.M()<<" correction "<<calc_wtag_corr(wjet.eta(), wjet.pt())<<" final mass "<<mass <<endl;
    //if(mass < 65 || mass > 95) continue;
    if(mass < 65 || mass > 105) continue;
    tmp_hyp.set_wHad(wjet.v4());
    for(auto & neutrino : neutrinos){
      if(deltaR(wjet.v4(),neutrino+lep)>dRmin){
	tmp_hyp.set_wLep(neutrino+lep);
	recoHyps.push_back(tmp_hyp);
      }
    }
  }
  
  int size_recoHyps = recoHyps.size();
  //cout<<"number of hyps "<<size_recoHyps<<endl;
  //cout<<"================================="<<endl;
  event.set(hypothesis,move(recoHyps));
  return size_recoHyps > 0 ? true : false;
}


double BprimeReco::calc_wtag_corr(double eta, double pt){
  
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
        
  genCorr =  puppisd_corrGEN->Eval( pt );
  if( fabs(eta)  <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( pt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( pt );
  }
  
  totalWeight = genCorr * recoCorr;

  return totalWeight;
}
