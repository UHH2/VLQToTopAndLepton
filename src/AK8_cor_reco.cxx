#include "UHH2/VLQToTopAndLepton/include/AK8_cor_reco.h"

AK8_cor_reco::AK8_cor_reco(uhh2::Context & ctx, std::string topjetcollection, TopJetId toptagid, TopJetId wjetid){
  is_mc = ctx.get("dataset_type") == "MC";
  if(!is_mc) return;
  h_topjets = ctx.get_handle<std::vector<TopJet> >(topjetcollection);
  for(unsigned int i=0; i<correction_names.size(); ++i){
    std::string name = correction_names[i];
    topjet_handles.push_back(ctx.declare_event_output<std::vector<TopJet>>("subjets_"+name));
    topjet_handles.push_back(ctx.declare_event_output<std::vector<TopJet>>("topjets_"+name));

    unsigned int m = i*2;
    std::unique_ptr<GenericTopJetCorrector> topjetcorr_help;
    topjetcorr_help.reset(new GenericTopJetCorrector(ctx,JERFiles::Summer16_23Sep2016_V4_L123_AK8PFpuppi_MC, "topjets_"+name, correction_int[i]));
    corrector_topjet.push_back(move(topjetcorr_help));
    std::unique_ptr<GenericSubJetCorrector> subjetcorr_help;
    subjetcorr_help.reset(new GenericSubJetCorrector(ctx,JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, "subjets_"+name, correction_int[i]));
    corrector_subjet.push_back(move(subjetcorr_help));
    
    reco.push_back(uhh2::make_unique<BprimeReco>(ctx,"wtag_reco_"+name));
    reco.push_back(uhh2::make_unique<BprimeReco>(ctx,"toptag_reco_"+name));
    reco[m+1]->set_topjetRecoId(toptagid);
    reco[m+1]->set_topjetCollection(ctx,"subjets_"+name);
    reco[m]->set_topjetRecoId(wjetid);
    reco[m]->set_topjetCollection(ctx,"topjets_"+name);

    dis.push_back(uhh2::make_unique<BprimeDiscriminator>(ctx,BprimeDiscriminator::wTag,"wtag_reco_"+name,"wtag_dis_"+name));
    dis.push_back(uhh2::make_unique<BprimeDiscriminator>(ctx,BprimeDiscriminator::cmsTopTag,"toptag_reco_"+name,"toptag_dis_"+name));

    dis[m]->set_emptyHyp(true);
    dis[m+1]->set_emptyHyp(true);
    
  }
}

bool AK8_cor_reco::process(uhh2::Event & event){
  if(event.isRealData) return true;
  std::vector<TopJet> jets;
  int counter =-1;
  for(auto & jetcoll : topjet_handles){
    event.set(jetcoll, jets);
    counter++;
    if(counter%2==0){
      //std::cout<<"copying toptag-jets. counter:"<<counter<<std::endl;
      for(unsigned int i=0; i<event.topjets->size();++i){
	const TopJet j = event.topjets->at(i);
	event.get(jetcoll).push_back(j);
      }
    }
    else{
      const auto wtagjets = &event.get(h_topjets);
      for(unsigned int i=0; i<wtagjets->size();++i){
	const TopJet j = wtagjets->at(i);
	event.get(jetcoll).push_back(j);
      }
    }
  }
  
  for(auto &c : corrector_topjet){
    c->process(event);
  }
  for(auto &c : corrector_subjet){
    c->process(event);
  }

  counter = -1;
  for(auto &r : reco){
    counter++;
    if(counter%2==0)
      r->hadronicW(event,2);
    else
      r->TopJetReco(event,2);
  }
  
  for(auto &d : dis){
    d->process(event);
  }
  
  return true;
}
