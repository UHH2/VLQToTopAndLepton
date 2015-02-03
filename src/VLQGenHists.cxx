//framework includes
#include "UHH2/VLQToTopAndLepton/include/VLQGenHists.h"


//#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/core/include/Event.h"

//general includes
#include <iostream>

using namespace std;
using namespace uhh2;


GenParticleHists VLQGenHists::histoBooker(string HistName, double minMass, double maxMass){
  GenParticleHists hists;
  
  hists.number = book<TH1F>((HistName+string("_number")).c_str(),(string("Number of ")+HistName).c_str(), 5, -0.5, 4.5);
  hists.decay_mom = book<TH1F>((HistName+string("_decay_mom")).c_str(), (HistName+string(" decay modes mom")).c_str(), 61, -30.5, 30.5);
  hists.decay_daughter = book<TH1F>((HistName+string("_decay_daughter")).c_str(), (HistName+string(" decay modes daughter")).c_str(), 61, -30.5, 30.5);
				    
  string suffix ="";
  string axisSuffix ="";

  for(int i = 0; i< 3; ++i){
    singleHists single; 
    if(i>0){ suffix = string("_")+to_string(i); axisSuffix = to_string(i);}
    single.pt  = book<TH1F>((HistName+string("_pt")+suffix).c_str(), (string("p_{T}^{")+HistName+string(" ")+axisSuffix+string("} [GeV/c]")).c_str(), 200, 0, 1500);
    single.eta = book<TH1F>((HistName+string("_eta")+suffix).c_str(),(string("#eta_{")+HistName+string(" ")+axisSuffix+string("}")).c_str(), 100, -4, 4);
    single.phi = book<TH1F>((HistName+string("_phi")+suffix).c_str(),(string("#phi_{")+HistName+string(" ")+axisSuffix+string("}")).c_str(), 100, -3.2, 3.2);
    single.mass = book<TH1F>((HistName+string("_mass")+suffix).c_str(),(string("mass_{")+HistName+string(" ")+axisSuffix+string("}")).c_str(), 100, minMass, maxMass);
    single.charge = book<TH1F>((HistName+string("_charge")+suffix).c_str(),(string("charge_{")+HistName+string(" ")+axisSuffix+string("}")).c_str(), 100, -1, 1);
    
    hists.stdHists.push_back(single);
    
  }
  
  return hists;
}

void VLQGenHists::histoFiller(vector<GenParticle> & particles, int partNumber, double weight){
  if(partNumber==-1) assert(0==1);

  sort_by_pt(particles);  
  m_Hists.at(partNumber).number->Fill(particles.size());

  int count = 0;
  for(auto part : particles){
    count++;
    for(int it = 0; it<3; ++it){
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).pt->Fill(part.pt(),weight);
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).eta->Fill(part.eta(),weight);
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).phi->Fill(part.phi(),weight);
      if(count==it || it==0 ) m_Hists.at(partNumber).stdHists.at(it).mass->Fill(part.v4().M(),weight);
      //if(count==it || it ==0) m_Hists.at(partNumber).stdHists.at(it).mass->Fill(sqrt(part.energy()*part.energy()-part.pt()*part.pt()),weight);
      if(count==it || it ==0) m_Hists.at(partNumber).stdHists.at(it).charge->Fill(part.charge(),weight);
    }
  }
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

int VLQGenHists::positionHelper(string Name){

  for(unsigned int i = 0; i < PartNames.size(); ++i ){
    if( PartNames.at(i).compare(Name)==0) return i; 
  }

  return -1;
}

VLQGenHists::VLQGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  
  vector<string> names {"VLQ","Higgs","Z","W","Top","B","Muon","Electron"};



  for (unsigned int i = 0; i< names.size(); ++i) {
    PartNames.push_back(names.at(i));
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

}

 


void VLQGenHists::fill(const Event & event){
  // fill the histograms. Don't forget to always use the weight when filling:
  //     double weight = event.weight;
  double weight = event.weight;
  
  // for (unsigned int i=0; i<event.get_current_triggernames().size();i++)
  //  cout<< event.get_current_triggernames()[i]<<endl;

     
  const vector<GenParticle> * genparticles = event.genparticles;
  
  vector<GenParticle> vlq;
  
  vector<GenParticle> higgs;
  vector<GenParticle> wboson;
  vector<GenParticle> zboson;
  
  vector<GenParticle> top;
  vector<GenParticle> bquarks;
  vector<GenParticle> electrons;
  vector<GenParticle> muons;
  
  
  for(auto igenp : *genparticles){ 
    
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



    if(abs(daughter1_pdgId)==6000008 ||abs(daughter1_pdgId)==6000007 || abs(daughter2_pdgId)==6000008 ||abs(daughter2_pdgId)==6000007){
      VLQ_mother->Fill(igenp.pdgId(),weight);
    }


    //put all the particles in vectors
    if (abs(igenp.pdgId()) == 6000008 || abs(igenp.pdgId()) == 6000007){
      vlq.push_back(igenp);
      decayFiller(weight,positionHelper("VLQ") , 0,0,daughter1_pdgId,daughter2_pdgId);
      VLQ_mother1_mother2->Fill(mother1_pdgId,mother2_pdgId,weight);
    }
    else if(abs(igenp.pdgId()) == 5 ){
      bquarks.push_back(igenp);
      decayFiller(weight,positionHelper("B") , 0,0,daughter1_pdgId,daughter2_pdgId);
    }
    else if(abs(igenp.pdgId()) == 6 ){
      top.push_back(igenp);     
      decayFiller(weight,positionHelper("Top") , 0,0,daughter1_pdgId,daughter2_pdgId);
    }
    else if (abs(igenp.pdgId()) == 11 ){
      electrons.push_back(igenp);
      decayFiller(weight,positionHelper("Electron") , 0,0,daughter1_pdgId,daughter2_pdgId);
    }
    else if (abs(igenp.pdgId()) == 13 ){
      muons.push_back(igenp);
      decayFiller(weight,positionHelper("Muon") , 0,0,daughter1_pdgId,daughter2_pdgId);
    }
    else if(abs(igenp.pdgId()) == 23 ){
      zboson.push_back(igenp);
      decayFiller(weight,positionHelper("Z") , 0,0,daughter1_pdgId,daughter2_pdgId);
    }
    else if(abs(igenp.pdgId()) == 24 ){
      wboson.push_back(igenp);
      decayFiller(weight,positionHelper("W") , 0,0,daughter1_pdgId,daughter2_pdgId);
    }
    else if(abs(igenp.pdgId()) == 25 ){
      higgs.push_back(igenp);
      decayFiller(weight,positionHelper("Higgs") , 0,0,daughter1_pdgId,daughter2_pdgId);
    }
    
    //if(abs(igenp.pdgId())==25 ||abs(igenp.pdgId())==23)
    //cout<<"pdgId GenParticle: "<<igenp.pdgId() <<" mom1: "<<  mother1_pdgId<<" mom2: "<< mother2_pdgId <<" daughter1: "<<  daughter1_pdgId<<" daughter2: "<< daughter2_pdgId<<endl;
  
    //Fill decay histograms

    
    if(abs(mother1_pdgId)== 6000007 || abs(mother1_pdgId)== 6000008 || abs(mother2_pdgId)== 6000008 || abs(mother2_pdgId)== 6000007 ){
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

  histoFiller(vlq, positionHelper("VLQ"), weight);
  histoFiller(higgs, positionHelper("Higgs"), weight);
  histoFiller(wboson, positionHelper("W"), weight);
  histoFiller(zboson, positionHelper("Z"), weight);
  histoFiller(top, positionHelper("Top"), weight);
  histoFiller(bquarks, positionHelper("B"), weight);
  histoFiller(electrons, positionHelper("Electron"), weight);
  histoFiller(muons, positionHelper("Muon"), weight);

 
}


VLQGenHists::~VLQGenHists(){}
