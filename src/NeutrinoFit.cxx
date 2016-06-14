#include "UHH2/VLQToTopAndLepton/include/NeutrinoFit.h"

using namespace std;
using namespace uhh2;

NeutrinoFit::NeutrinoFit(bool fit){
  do_fit = fit;
  positiv = new TMinuit(5);
}
NeutrinoFit::~NeutrinoFit(){}

double NeutrinoFit::DeltaPolarNeutrino(double PhiN, double metPx, double metPy, double PhiLep, double PtLep){
  double PyN;
  double PxN;
  const double mass_w = 80.399;
  double PI = boost::math::constants::pi<double>();

  double num = 10.e-7;

  if(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN )) < num){
    PyN = 0.5*mass_w*mass_w* sin(PhiN)/(PtLep*num);
    PxN = 0.5*mass_w*mass_w* cos(PhiN)/(PtLep*num);
  }
  else{
    PyN = 0.5*mass_w*mass_w* sin(PhiN)/(PtLep*(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN ))));
    PxN = 0.5*mass_w*mass_w* cos(PhiN)/(PtLep*(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN ))));
  }

  return pow(PxN-metPx,2)+pow(PyN-metPy,2);

}

void NeutrinoFit::polarminuitfunc(int& nDim, double* gout, double& result, double par[], int flg){
  result = DeltaPolarNeutrino(par[0],par[1],par[2],par[3],par[4]);
}

std::vector<LorentzVector> NeutrinoFit::NeutrinoFitPolar(const LorentzVector lepton, const LorentzVector met){
  double PI = boost::math::constants::pi<double>();
  TVector3 lepton_pT = toVector(lepton);
  lepton_pT.SetZ(0);

  TVector3 neutrino_pT = toVector(met);
  neutrino_pT.SetZ(0);

  const double mass_w = 80.399;

  double min = -2*PI;
  double max = 2*PI;
  double start = met.phi();
  double step = 10.e-5;

  double mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
  double A = - (lepton_pT * lepton_pT);
  double B = mu * lepton.pz();
  double C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);

  double discriminant = B * B - A * C;

  std::vector<LorentzVector> solutions;

  if (0 >= discriminant && do_fit){

    double resultPhi = 0;
    double error = 0;
    int ierflg;
    double* arg = new double[1];

    positiv->SetPrintLevel(-1); // -1 quiet, 0 normal, 1 verbose; Preset 0

    positiv->SetFCN(polarminuitfunc);

    positiv->DefineParameter(0,"PhiN",start, step,  min, max);
    positiv->DefineParameter(1,"metPx",met.px(),0,0,0);
    positiv->DefineParameter(2,"metPy",met.py(),0,0,0);
    positiv->DefineParameter(3,"PhiLep",lepton.phi(),0,0,0);
    positiv->DefineParameter(4,"PtLep",lepton.pt(),0,0,0);

    positiv->FixParameter(1);
    positiv->FixParameter(2);
    positiv->FixParameter(3);
    positiv->FixParameter(4);

    positiv->SetMaxIterations(500);

    arg[0]= 2;
    positiv->mnexcm("SET STR",arg,1,ierflg);

    positiv->Migrad();

    positiv->GetParameter(0,resultPhi,error);

    delete[] arg;

    if(resultPhi != resultPhi){
	std::cerr << "neutrino phi is NAN " << std::endl;
    }
    if(resultPhi > PI) resultPhi = resultPhi-2*PI;
    if(resultPhi < PI) resultPhi = resultPhi+2*PI;
    double PyN;
    double PxN;

    double num = 10.e-7;

    if(1-cos(deltaPhiAbs(lepton.phi(), resultPhi)) < num){
      PyN = 0.5*mass_w*mass_w* sin(resultPhi)/(lepton.pt()*num);
      PxN = 0.5*mass_w*mass_w* cos(resultPhi)/(lepton.pt()*num);
    }
    else{
      PyN = 0.5*mass_w*mass_w* sin(resultPhi)/(lepton.pt()*(1-cos(deltaPhiAbs(lepton.phi(), resultPhi))));
      PxN = 0.5*mass_w*mass_w* cos(resultPhi)/(lepton.pt()*(1-cos(deltaPhiAbs(lepton.phi(), resultPhi))));
    }

    LorentzVectorXYZE neutrino_result(0,0,0,0);
    neutrino_result.SetPx(PxN);
    neutrino_result.SetPy(PyN);

    double pzfit =  lepton.pz()*neutrino_result.pt()/lepton.pt();

    LorentzVectorXYZE solution (0,0,0,0);
    solution.SetPx(PxN);
    solution.SetPy(PyN);
    solution.SetPz(pzfit);
    solution.SetE(toVector(solution).Mag());

    solutions.push_back(toPtEtaPhi(solution));
  }
  else{
    if(discriminant < 0)
      discriminant = 0; 

    discriminant = sqrt(discriminant);
    
    LorentzVectorXYZE solution (0,0,0,0);
    solution.SetPx(met.Px());
    solution.SetPy(met.Py());
    solution.SetPz((-B - discriminant) / A);
    solution.SetE(toVector(solution).Mag());

    solutions.push_back(toPtEtaPhi(solution));
    if(discriminant==0)
      return solutions;

    LorentzVectorXYZE solution2 (0,0,0,0);
    solution2.SetPx(met.Px());
    solution2.SetPy(met.Py());
    solution2.SetPz((-B + discriminant) / A);
    solution2.SetE(toVector(solution2).Mag());

    solutions.push_back(toPtEtaPhi(solution2));
  }
  return solutions;
}

double NeutrinoFit::deltaPhiAbs(double x1, double x2)
{
  // x1 & x2 are two Phi expected in the range [-PI,PI]
  double deltaphi = fabs(x1 - x2);
  if(deltaphi > M_PI) deltaphi = 2*M_PI - deltaphi;
  return deltaphi;
}
