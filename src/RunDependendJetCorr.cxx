#include "UHH2/VLQToTopAndLepton/include/RunDependendJetCorr.h"



RunDependendJetCorr::RunDependendJetCorr(uhh2::Context & ctx, std::string topjetcollection, std::string subjetcollection, int direction){
  is_mc = ctx.get("dataset_type") == "MC";
  if(subjetcollection.empty())subjetcollbool=false; 
  if(is_mc){
    jet_corrector_MC.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK8PFpuppi_MC,topjetcollection, direction));
    //subjet_corrector_MC.reset(new SubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, direction));
    if(subjetcollection.empty())subjet_corrector_MC.reset(new SubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, direction));
    else subjet_corrector_MC_coll.reset(new GenericSubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC,subjetcollection, direction));
  }
  else{
    jet_corrector_BCD.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFpuppi_DATA,topjetcollection));
    jet_corrector_EFearly.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFpuppi_DATA,topjetcollection));
    jet_corrector_FlateG.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFpuppi_DATA,topjetcollection));
    jet_corrector_H.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFpuppi_DATA,topjetcollection));
    
    subjet_corrector_BCD.reset(new SubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA));
    subjet_corrector_EFearly.reset(new SubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA));
    subjet_corrector_FlateG.reset(new SubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA));
    subjet_corrector_H.reset(new SubJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA));
  }



  
}

bool RunDependendJetCorr::process(uhh2::Event & event){
  if(is_mc){
    jet_corrector_MC->process(event);
    if(!subjetcollbool) subjet_corrector_MC->process(event);
    else subjet_corrector_MC_coll->process(event);
  }
  else{
    if(event.run <= runnr_BCD) {
      jet_corrector_BCD->process(event);
      jet_corrector_BCD->process(event);
    }
    else if(event.run < runnr_EFearly){
      jet_corrector_EFearly->process(event); //< is correct, not <=
      subjet_corrector_EFearly->process(event); //< is correct, not <=
	    
    }
    else if(event.run <= runnr_FlateG){
      jet_corrector_FlateG->process(event);
      subjet_corrector_FlateG->process(event);
    }
    else if(event.run > runnr_FlateG){
      jet_corrector_H->process(event);
      subjet_corrector_H->process(event);
    }
    else throw std::runtime_error("RunDependendJetCorr.cxx: run number not covered by if-statements in process-routine.");
  }
  return true;
}


namespace JERFiles {

  const std::vector<std::string> Summer16_23Sep2016_V4_L123_AK8PFpuppi_MC = {
    "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt",
  };
  
  const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK8PFpuppi_DATA = {
    "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt",
  };
  
  const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK8PFpuppi_DATA = {
    "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK8PFchs.txt",
  };
  
  const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK8PFpuppi_DATA = {
    "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK8PFchs.txt",
  };
  
  const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK8PFpuppi_DATA = {
    "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK8PFchs.txt",
    "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK8PFchs.txt",
  };
}
