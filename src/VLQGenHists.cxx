//framework includes
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"


//#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/core/include/Event.h"

//general includes
#include <iostream>

using namespace std;
using namespace uhh2;


GenParticleHists VLQGenHists::histoBooker(const string& HistName, double minMass, double maxMass){
  GenParticleHists hists;
  
  hists.number = book<TH1F>(HistName+"_number","Number of "+HistName, 5, -0.5, 4.5);
  hists.decay_mom = book<TH1F>(HistName+"_decay_mom", HistName+" decay modes mom", 61, -30.5, 30.5);
  hists.decay_daughter = book<TH1F>(HistName+"_decay_daughter", HistName+" decay modes daughter", 61, -30.5, 30.5);
  hists.deltaRmin = book<TH1F>(HistName+"_deltaRmin","DeltaRmin of "+HistName, 40, 0, 8);
  hists.nextParticle = book<TH1F>(HistName+"_nextParticle","Nearest particle of "+HistName, 61, -30.5,30.5);
  hists.deltaR_daughters = book<TH1F>(HistName+"_deltaR_daughters","DeltaR daughters "+HistName, 100, 0,10);
  
		    
  string suffix ="";
  string axisSuffix ="";

  for(int i = 0; i< 3; ++i){
    singleHists single; 
    if(i>0){ suffix = "_"+to_string(i); axisSuffix = to_string(i);}

    single.pt  = book<TH1F>(HistName+"_pt"+suffix, "p_{T}^{"+HistName+" "+axisSuffix+"} [GeV/c]", 200, 0, 1500);
    single.eta = book<TH1F>(HistName+"_eta"+suffix,"#eta_{"+HistName+" "+axisSuffix+"}", 100, -5, 5);
    single.phi = book<TH1F>(HistName+"_phi"+suffix,"#phi_{"+HistName+" "+axisSuffix+"}", 100, -3.2, 3.2);
    single.mass = book<TH1F>(HistName+"_mass"+suffix,"mass_{"+HistName+" "+axisSuffix+"}", 100, minMass, maxMass);
    single.charge = book<TH1F>(HistName+"_charge"+suffix,"charge_{"+HistName+" "+axisSuffix+"}", 100, -5, 5);
    single.pt_eta = book<TH2F>(HistName+"_pT_eta"+suffix,"pT & eta "+HistName+" "+axisSuffix, 150, 0, 800, 100, -4, 4);

    hists.stdHists.push_back(single);
    
  }
  
  return hists;
}

void VLQGenHists::histoFiller(vector<GenParticle> & particles,  int partNumber, double weight){
  if(partNumber==-1) assert(0==1);

  sort_by_pt(particles);  
  m_Hists.at(partNumber).number->Fill(particles.size(),weight);

  int count = 0;
  for(auto part : particles){
    count++;
    for(int it = 0; it<3; ++it){
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).pt->Fill(part.pt(),weight);
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).eta->Fill(part.eta(),weight);
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).phi->Fill(part.phi(),weight);
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).mass->Fill(part.v4().M(),weight);
      if(count==it || it ==0) m_Hists.at(partNumber).stdHists.at(it).charge->Fill(part.charge(),weight);
      //if(partNumber==6)cout<<" num "<<it<<" "<<part.charge()<<endl;
      if(count==it || it ==0) m_Hists.at(partNumber).stdHists.at(it).pt_eta->Fill(part.pt(),part.eta(),weight);

    }
  }
}
void VLQGenHists::fill_nearest(int position, double weight, double deltaRmin, double pdgId){
  m_Hists.at(position).deltaRmin->Fill(deltaRmin,weight);
  m_Hists.at(position).nextParticle->Fill(pdgId,weight);
}

void VLQGenHists::fillDaughterDistance(int position, double deltaR, double weight){
  m_Hists.at(position).deltaR_daughters->Fill(deltaR,weight);
}

void VLQGenHists::decayFiller(double weight, int partNumber, int mother1, int mother2, int daughter1, int daughter2){
 
  if(mother1!=0)
    m_Hists.at(partNumber).decay_mom->Fill(mother1,weight);
  if(mother2!=0)
    m_Hists.at(partNumber).decay_mom->Fill(mother2,weight);
  if(daughter1!=0)
    m_Hists.at(partNumber).decay_daughter->Fill(daughter1,weight);
  if(daughter2!=0)
    m_Hists.at(partNumber).decay_daughter->Fill(daughter2,weight);
}

int VLQGenHists::positionHelper(const string& Name){
  for(unsigned int i = 0; i < PartNames.size(); ++i ){
    if( PartNames.at(i).compare(Name)==0) return i; 
  }
  return -1;
}



VLQGenHists::VLQGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  //make sure that these two vectors corresspond, vlq are not considered
  vector<string> names {"VLQ","Higgs","Z","W","Top","B","Muon","muNeutrino","Electron","eleNeutrino","forward"};
  vector<int> pdgIds {25,23,24,6,5,13,14,11,12};  

   for (unsigned int i = 0; i< names.size(); ++i) {
    PartNames.push_back(names.at(i));
    if(i<pdgIds.size())PartPdgId.push_back(pdgIds.at(i));
    double minMass = 0;
    double maxMass = 300;
    if(i==0) maxMass = 1600;
    if(i>4) maxMass = 10;
    m_Hists.push_back(histoBooker(PartNames[i],minMass,maxMass));
    maxMass =300;
  }

  //no mothers
  particles_noMother     = book<TH1F>("particle_noMother"   , "Particles with no mother", 61, -30.5, 30.5);
  particles_noMother_pT  = book<TH1F>("particles_noMother_pT" , "p_{T} [GeV/c]", 200, 0, 1000);
  particles_noMother_eta = book<TH1F>("particles_noMother_eta", "#eta ", 40, -2.5, 2.5);
  particles_noMother_phi = book<TH1F>("particles_noMother_phi", "#phi", 64, -3.2, 3.2);
  //no daughters or mothers
  particles_noMotherNoDaughter     = book<TH1F>("particle_noMotherNoDaughter"   , "Particles with no mother or daughter", 61, -30.5, 30.5);
  particles_noMotherNoDaughter_pT  = book<TH1F>("particles_noMotherNoDaughter_pT" , "p_{T} [GeV/c]", 200, 0, 1000);
  particles_noMotherNoDaughter_eta = book<TH1F>("particles_noMotherNoDaughter_eta", "#eta ", 40, -2.5, 2.5);
  particles_noMotherNoDaughter_phi = book<TH1F>("particles_noMotherNoDaughter_phi", "#phi", 64, -3.2, 3.2);

  VLQ_mother = book<TH1F>("VLQ_mothers"   , "VLQ mothers", 61, -30.5, 30.5);
  VLQ_mother1_mother2= book<TH2F>("VLQ_mother1_mother2"   , "VLQ mothers", 61, -30.5, 30.5,61, -30.5, 30.5);
  VLQ_E_daughters = book<TH2F>("VLQ_E_daughters"   , "VLQ E daughters", 100, 0, 3000, 100, 0, 3000);
  VLQ_Eratio_daughters = book<TH1F>("VLQ_Eratio_daughters"   , "VLQ E_{ratio} daughters", 100, 0, 5);
  VLQ_pt_daughters = book<TH2F>("VLQ_pt_daughters"   , "VLQ pt daughters", 100, 0, 1000, 100, 0, 1000);
  VLQ_ptratio_daughters = book<TH1F>("VLQ_ptratio_daughters"   , "VLQ pt_{ratio} daughters", 100, 0, 5);

  deltaR_w   = book<TH1F>("deltaR_w "   , "#Delta R(W_{1},W_{2})", 100, 0, 8);
  deltaPhi_w = book<TH1F>("deltaPhi_w"  , "#Delta #phi(W_{1},W_{2})", 100, 0, 4);

   deltaPhi_t_w = book<TH1F>("deltaPhi_t_w","#delta #phi (t,W)", 100, 0, 8);
   deltaEta_t_w = book<TH1F>("deltaEta_t_w","#delta #eta (t,W)", 100, 0, 8);
   deltaPhi_deltaEta_t_w = book<TH2F>("deltaPhi_deltaEta_t_w","", 100, 0, 8, 100,0,8);  
}


void VLQGenHists::fill(const Event & event){

  if(event.isRealData) return;
  // fill the histograms. Don't forget to always use the weight when filling:
  //     double weight = event.weight;
  double weight = event.weight;
  const vector<GenParticle> * genparticles = event.genparticles;
  vector<GenParticle> vlq;
  vector<GenParticle> forwardJet;
  vector<vector<GenParticle>> particleStore;
  
  for(unsigned int i =0; i<PartPdgId.size();++i){
    vector<GenParticle> oneType;
    particleStore.push_back(oneType);
  }
  
  for(auto & igenp : *genparticles){  
    const GenParticle* mother1 = igenp.mother(genparticles,1);
    const GenParticle* mother2 = igenp.mother(genparticles,2);
    
    int mother1_pdgId = 0; 
    if(mother1)mother1_pdgId = mother1->pdgId();
    int mother2_pdgId = 0;
    if(mother2)mother2->pdgId();
    
    const GenParticle* daughter1 = igenp.daughter(genparticles,1);
    const GenParticle* daughter2 = igenp.daughter(genparticles,2);
    
    int daughter1_pdgId = 0; 
    if(daughter1)daughter1_pdgId = daughter1->pdgId();
    int daughter2_pdgId = 0;
    if(daughter2) daughter2_pdgId = daughter2->pdgId();
    
    auto nearGen =  closestParticle_noNeutrino(igenp, *genparticles);

    for(unsigned int i =0; i<PartPdgId.size();++i){
      if(abs(igenp.pdgId())<5 ||abs(igenp.pdgId())>25)break;
      if(PartPdgId.at(i) == abs(igenp.pdgId())){
	particleStore.at(i).push_back(igenp);
        decayFiller(weight,i+1,0,0,daughter1_pdgId,daughter2_pdgId);
        fill_nearest(i+1,weight,uhh2::deltaR(*nearGen,igenp),nearGen->pdgId());
	if(daughter1 && daughter2 )fillDaughterDistance(i+1, uhh2::deltaR(*daughter1,*daughter2),weight);
	break;
      }
    }
    
    if(abs(daughter1_pdgId)>6000000 || abs(daughter2_pdgId)>6000000){
      VLQ_mother->Fill(igenp.pdgId(),weight);
    }
    //Energy & pt ratios of daughters of VLQ
    if(abs(igenp.pdgId())>6000000){
      //cout<<"Found B'"<<endl;
      //cout<<abs(daughter1_pdgId)<<" "<<abs(daughter2_pdgId)<<endl;
      if(abs(daughter1_pdgId) > 0 && abs(daughter2_pdgId) > 0){
	//cout<<abs(daughter1_pdgId)<<" "<<abs(daughter2_pdgId)<<endl;
	//cout<<"Gen t,W delta phi "<<deltaPhi(daughter1->v4(),daughter2->v4())<<" delta eta "<<abs(daughter2->v4().eta()-daughter1->v4().eta())<<endl; 
	deltaPhi_t_w->Fill(deltaPhi(daughter1->v4(),daughter2->v4()),weight);
	deltaEta_t_w->Fill(abs(daughter2->v4().eta()-daughter1->v4().eta()),weight);
	deltaPhi_deltaEta_t_w->Fill(deltaPhi(daughter1->v4(),daughter2->v4()),abs(daughter2->v4().eta()-daughter1->v4().eta()),weight);
	if(abs(daughter1_pdgId)==6){
	  VLQ_E_daughters->Fill(daughter1->energy(),daughter2->energy(),weight);
	  VLQ_Eratio_daughters->Fill(daughter1->energy()/daughter2->energy(),weight);
	  VLQ_pt_daughters->Fill(daughter1->pt(),daughter2->pt(),weight);
	  VLQ_ptratio_daughters->Fill(daughter1->pt()/daughter2->pt(),weight);
	}
	else if(abs(daughter1_pdgId)==24){
	  VLQ_E_daughters->Fill(daughter2->energy(),daughter1->energy(),weight);
	  VLQ_Eratio_daughters->Fill(daughter2->energy()/daughter1->energy(),weight);
	  VLQ_pt_daughters->Fill(daughter2->pt(),daughter1->pt(),weight);
	  VLQ_ptratio_daughters->Fill(daughter2->pt()/daughter1->pt(),weight);
	}
      }
    }
    //forwardJet
    if(abs(daughter1_pdgId) < 5 && abs(daughter1_pdgId) < 0 && abs(daughter2_pdgId) > 6000000)forwardJet.push_back(*daughter1);
    if(abs(daughter2_pdgId) < 5 && abs(daughter2_pdgId) < 0 && abs(daughter1_pdgId) > 6000000)forwardJet.push_back(*daughter2);
    
    if(abs(igenp.pdgId()) <5){
      if((abs(mother1_pdgId) <5 || abs(mother1_pdgId) ==21)  && abs(mother2_pdgId) == 0 && abs(daughter1_pdgId)==0 &&  abs(daughter2_pdgId)==0)
	forwardJet.push_back(igenp);
      else if(abs(mother1_pdgId) ==0 && (abs(mother2_pdgId) <5 || abs(mother2_pdgId) ==21) && abs(daughter1_pdgId)==0 &&  abs(daughter2_pdgId)==0)
	forwardJet.push_back(igenp);
    }
    /*
    if(abs(igenp.pdgId()) <5 && (abs(mother1_pdgId) == 24 || abs(mother2_pdgId) == 24)){
      auto help1_mother1 = mother1->mother(genparticles,1)->pdgId();
      auto help1_mother2 = mother1->mother(genparticles,2)->pdgId();
      auto help2_mother1 = mother2->mother(genparticles,1)->pdgId();
      auto help2_mother2 = mother2->mother(genparticles,2)->pdgId();
      if(abs(help1_mother1) != 6 && abs(help1_mother2) != 6 && abs(help2_mother1) != 6 && abs(help2_mother2) != 6 )
	forwardJet.push_back(igenp);
    }
    */
 
    //put all the particles in vectors
    if (abs(igenp.pdgId())> 6000000){
      vlq.push_back(igenp);
      decayFiller(weight,positionHelper("VLQ"),0,0,daughter1_pdgId,daughter2_pdgId);
      VLQ_mother1_mother2->Fill(mother1_pdgId,mother2_pdgId,weight);
      if(daughter1 && daughter2 )fillDaughterDistance(0, uhh2::deltaR(*daughter1,*daughter2),weight);
      fill_nearest(0,weight, uhh2::deltaR(*nearGen,igenp),nearGen->pdgId());  
    }
    
    //if(abs(igenp.pdgId())==25 ||abs(igenp.pdgId())==23)
    //cout<<"pdgId GenParticle: "<<igenp.pdgId() << " eta "<<igenp.eta() <<" mom1: "<<  mother1_pdgId<<" mom2: "<< mother2_pdgId <<" daughter1: "<<  daughter1_pdgId<<" daughter2: "<< daughter2_pdgId<<endl;
    
    //Fill decay histograms  
    if(abs(mother1_pdgId)> 6000000 || abs(mother2_pdgId)> 6000000){
      decayFiller(weight,positionHelper("VLQ") , mother1_pdgId, mother2_pdgId,0,0);    
    }
    if(abs(mother1_pdgId)== 5 || abs(mother2_pdgId)== 5)
      decayFiller(weight,positionHelper("B") , mother1_pdgId, mother2_pdgId,0,0);
    if(abs(mother1_pdgId)== 6 || abs(mother2_pdgId)== 6)
      decayFiller(weight,positionHelper("Top") , mother1_pdgId, mother2_pdgId,0,0); 
    if(abs(mother1_pdgId)== 11 || abs(mother2_pdgId)== 11)
      decayFiller(weight,positionHelper("Electron") , mother1_pdgId, mother2_pdgId,0,0);
    if(abs(mother1_pdgId)== 13 || abs(mother2_pdgId)== 13)
      decayFiller(weight,positionHelper("Muon") , mother1_pdgId, mother2_pdgId,0,0);
    if(abs(mother1_pdgId)== 23 || abs(mother2_pdgId)== 23)
      decayFiller(weight,positionHelper("Z") , mother1_pdgId, mother2_pdgId,0,0); 
    if(abs(mother1_pdgId)== 24 || abs(mother2_pdgId)== 24)
      decayFiller(weight,positionHelper("W") , mother1_pdgId, mother2_pdgId,0,0); 
    if(abs(mother1_pdgId)== 25 || abs(mother2_pdgId)== 25)
      decayFiller(weight,positionHelper("Higgs") , mother1_pdgId, mother2_pdgId,0,0); 
    if(abs(mother1_pdgId)== 12 || abs(mother2_pdgId)== 12)
      decayFiller(weight,positionHelper("eleNeutrino") , mother1_pdgId, mother2_pdgId,0,0);
    if(abs(mother1_pdgId)== 14 || abs(mother2_pdgId)== 14)
      decayFiller(weight,positionHelper("muNeutrino") , mother1_pdgId, mother2_pdgId,0,0);
    if(abs(mother1_pdgId)==5  || abs(mother2_pdgId)==5 )
      decayFiller(weight,positionHelper("B") , mother1_pdgId, mother2_pdgId,0,0);


    if(abs(mother1_pdgId)== 0 && abs(mother2_pdgId)== 0){
      particles_noMother->Fill(igenp.pdgId());           
      particles_noMother_pT->Fill(igenp.pt());  
      particles_noMother_eta->Fill(igenp.eta());  
      particles_noMother_phi->Fill(igenp.phi());   
    }
    if(abs(mother1_pdgId)== 0 && abs(mother2_pdgId)== 0 &&  abs(daughter2_pdgId)== 0 &&  abs(daughter1_pdgId)== 0){
      particles_noMotherNoDaughter->Fill(igenp.pdgId());           
      particles_noMotherNoDaughter_pT->Fill(igenp.pt());  
      particles_noMotherNoDaughter_eta->Fill(igenp.eta());  
      particles_noMotherNoDaughter_phi->Fill(igenp.phi()); 
      
    } 
  }
  /*
  cout<<forwardJet.size()<<endl;
  cout<<"--------------------------------------"<<endl;
  cout<<"--------------------------------------"<<endl;
  */

  histoFiller(vlq, positionHelper("VLQ"), weight);
  
  for(unsigned int i =0; i<PartPdgId.size();++i)
    histoFiller(particleStore.at(i),i+1,weight); 

  //cout<<positionHelper("forward")<<endl;
  histoFiller(forwardJet,positionHelper("forward"),weight);

  if(particleStore.at(positionHelper("W")-1).size()==2){
    deltaR_w->Fill(deltaR(particleStore.at(positionHelper("W")-1).at(0),particleStore.at(positionHelper("W")-1).at(1)),weight);
    deltaPhi_w->Fill(deltaPhi(particleStore.at(positionHelper("W")-1).at(0),particleStore.at(positionHelper("W")-1).at(1)),weight);
  }

}


VLQGenHists::~VLQGenHists(){}


const GenParticle* VLQGenHists::closestParticle_noNeutrino(const Particle  & p, const std::vector<GenParticle> & particles){
    double deltarmin = std::numeric_limits<double>::infinity();
    const GenParticle* next=0;
    for(unsigned int i=0; i<particles.size(); ++i) {
        const GenParticle & pi = particles[i];
        double dr = uhh2::deltaR(pi, p);
        if(dr < deltarmin && &pi != &p && pi.pdgId() != 12 && pi.pdgId() != 14 && pi.pdgId() != 16) {
            deltarmin = dr;
            next = &pi;
        }
    }
    return next;
}
