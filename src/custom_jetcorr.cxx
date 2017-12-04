#include "UHH2/VLQToTopAndLepton/include/custom_jetcorr.h"


custom_jetcorr::custom_jetcorr(uhh2::Context& ctx){
  JetId wide_softjet =  AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30, 5.0));
  int i =-1;
  std::vector<std::string> empty_string_vec = {};
  
  for(auto name : combination_names){
    i++;
    jet_handles.push_back(ctx.declare_event_output<std::vector<Jet>>("jet_"+name));
    met_handles.push_back(ctx.declare_event_output<MET>("met_"+name));
    
    jer_correctors.push_back(uhh2::make_unique<JetResolutionSmearer>(ctx,JERSmearing::SF_13TeV_2016,"jet_"+name,combination_ints[i][0]));
    jec_correctors.push_back(uhh2::make_unique<JetCorrector>(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, empty_string_vec,"jet_"+name,"met_"+name,combination_ints[i][1]));    
    jet_cleaners.push_back(uhh2::make_unique<JetCleaner>(ctx,wide_softjet,"jet_"+name));
  }
}

void custom_jetcorr::copy_jets(uhh2::Event & event){
  //copy ak4 jets
  std::vector<Jet> jets;
  for(auto & jetcoll : jet_handles)
    event.set(jetcoll, jets);
  for(unsigned int i=0; i<event.jets->size();++i){
    const Jet j = event.jets->at(i);
    for(auto & jetcoll : jet_handles)
      event.get(jetcoll).push_back(j);
  }
}

void custom_jetcorr::apply_corr(uhh2::Event & event){
  for(auto & jer_corr : jer_correctors)
    jer_corr->process(event);
  for(auto & jec_corr : jec_correctors)
    jec_corr->process(event);
}

void custom_jetcorr::clean(uhh2::Event & event){
  for(auto & jet_cleaner : jet_cleaners)
    jet_cleaner->process(event);
}

void custom_jetcorr::set_met(uhh2::Event & event){
  for(auto & handle : met_handles)
    event.set(handle, *event.met);
}
