////////////////////////////////////////////
// 
// RooParticle Test Application
// (Fixed Target) Rutherford Experiment
// 
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"

#include "iostream"
#include <math.h>
#include "unistd.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

class FixedTargetExperiment{
public:
  void setTemplate(event* eventT_){eventT=eventT_;}
  void makeEvent(particle* p); // if event already on,
  void offEvent(){eventOn=false;}
  void delEvent(){eventOn=false; delete currentEvent;}
  event* getEvent(){return currentEvent;}
protected:
private:
  event * eventT;
  Bool_t eventOn = false;
  event * currentEvent;
};

void FixedTargetExperiment::makeEvent(particle* p){
  eventOn = true;
  currentEvent = new event(*eventT);
  currentEvent->AddParticle(p);
}


// Float_t RandomFloat(Float_t a, Float_t b) {
//     Float_t random = ((Float_t) rand()) / (Float_t) RAND_MAX;
//     Float_t diff = b - a;
//     Float_t r = random * diff;
//     return a + r;
// }

TVector3 unitv(TVector3 v){
  TVector3 result = TVector3();
  result = v*(1/(v.Mag()));
  return result;
}

Double_t angleXD(TVector3 v){
  TVector3 unit = TVector3();
  unit = unitv(v);
  Double_t result = acos(unit[0]);
  return result;
}

int Experiment(Int_t nparticle, Int_t randomseed, TString outputFilename, Double_t derivingDegreeCriterion, Int_t derivingMinimumPoint, Double_t velocity, Double_t startingMinimumDistance, Int_t vperiod, Int_t autosavePeriod){
  event * eventT = new event();

  coulombForce * cp = new coulombForce();
  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  Double_t x1_[] = {0., 0., 0.};
  Double_t v1_[] = {0., 0., 0.};
  TVector3 x1 = TVector3(x1_);
  TVector3 v1 = TVector3(v1_);
  EMparticle * p1 = new EMparticle(206,79, x1, v1, true, false);

  eventT->AddParticle(p1);
  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  FixedTargetExperiment * FTE = new FixedTargetExperiment();
  FTE->setTemplate(eventT);

  TVector3 x2 = TVector3();
  TVector3 v2 = TVector3();

  Float_t imp_min = 0.;
  Float_t imp_max = +3000.;
  Float_t imp;

  EMparticle * p2;

  Double_t v2_[] = {velocity,0.,0.};
  v2 = TVector3(v2_);

  Double_t x2_[] = {-5000., 0., 0.,};

  Int_t i;
  Int_t j=0;

  TString fname = outputFilename;

  TFile * hfile = new TFile(fname,"RECREATE");
  TTree * tree = new TTree(TString::Format("TT_%d",j), "Rutherford Scattering Angle");

  std::cout<< "ROOT File save in " << fname << std::endl;

  Double_t   imp_tree;
  Int_t      NPOINT_tree;
  Double_t   fx_tree;
  Double_t   fy_tree;
  Double_t   sx_tree;
  Double_t   sy_tree;
  Double_t   vxF_tree;
  Double_t   vyF_tree;
  Double_t   SAngle_tree;
  Double_t   DCA_tree;
  Double_t   outtime_free;


  tree->Branch("imp", &imp_tree, "imp/D");
  tree->Branch("NPOINT", &NPOINT_tree, "Npoint/I");
  tree->Branch("fx", & fx_tree, "fx/D");
  tree->Branch("fy", & fy_tree, "fy/D");
  tree->Branch("sx", & sx_tree, "sx/D");
  tree->Branch("sy", & sy_tree, "sy/D");
  tree->Branch("vxF", & vxF_tree, "vxF/D");
  tree->Branch("vyF", & vyF_tree, "vyF/D");
  tree->Branch("SAngle", & SAngle_tree, "angle/D");
  tree->Branch("DCA", & DCA_tree, "DCA/D");

  TRandom3 * r1 = new TRandom3(randomseed);

  Double_t ImpWall = startingMinimumDistance * TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());

  for(i=0; i<nparticle ;i++ ){

    imp = r1->Uniform(imp_min, imp_max);

    x2_[1] = imp;
    if(imp>ImpWall){
      x2_[0] = -1 * imp / TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());
    }else{
      x2_[0] = -1 * startingMinimumDistance;
    }

    x2 = TVector3(x2_);
    p2 = new EMparticle(4,2, x2, v2, false, true);

    FTE->makeEvent(p2);

    FTE->getEvent()->AddInspector(new inspector(FTE->getEvent()->getParticle(0), p2, "DEG", derivingDegreeCriterion));
    FTE->getEvent()->AddInspector(new inspector(FTE->getEvent()->getParticle(0), p2, "CNT", derivingMinimumPoint));
    
    FTE->getEvent()->DeriveInspect(1,false);
    // FTE->getEvent()->DeriveDTN(1,10);


    imp_tree = imp;
    sx_tree = x2_[0];
    sy_tree = x2_[1];
    NPOINT_tree = p2->GetPath()->GetMaxNumber();
    fx_tree  = p2->GetPath()->GetLastX().operator[](0);
    fy_tree  = p2->GetPath()->GetLastX().operator[](1);
    vxF_tree = p2->GetPath()->GetLastV().operator[](0);
    vyF_tree = p2->GetPath()->GetLastV().operator[](1);
    SAngle_tree = angleXD(p2->GetPath()->GetLastV());
    p2->GetPath()->GetDCA(x1, DCA_tree, outtime_free);

    tree->Fill();

    if(i%autosavePeriod==0){
      tree->AutoSave();
    }

    if(i%vperiod==0){
      std::cout<<j<<"/"<<i<<" | IMP: "<<imp_tree<<" | xf_x: " << fx_tree << " | vf_x : "<<vxF_tree<<" | SAngle (DEG): "<<SAngle_tree * TMath::RadToDeg() <<" | POINT: "<<NPOINT_tree<<endl;
      
    }
    FTE->offEvent();
    FTE->delEvent();

  }

  tree->Print();
  tree->AutoSave();
  hfile->Close();
  delete hfile;

  return 0;
}

//===============================================================================
// ./test_rutherford_inspect
// Rutherford Experiment Simulation (Numerical)
// If there's less argv than parameters, it will use default value for absents
// -n   Number of MonteCarlo Simulation
// -r   Random Seed (UInt)
// -o   Output Filename (path)
// -d   Asymptote Degree Criterion
// -D   Minimum Incident Distance
// -p   Minumum Number of Derived Points
// -V   Incident Particle's Velocity (r.SOL)
// -v   Verbose Period
// -S   Autosave Period
// -h   Help
//===============================================================================

int main(int argc, char *argv[]){
  char c; // option

  //Default setting
  Int_t    nparticle               = 1000;
  Int_t    randomseed              = 65539;
  TString  outputFilename          = "Data/test_rutherford_inspect.root";
  Double_t derivingDegreeCriterion = 1;
  Int_t    derivingMinimumPoint    = 50000;
  Double_t startingMinimumDistance = 5000;
  Double_t velocity                = 0.05;
  Int_t    vperiod                 = 5;
  Int_t    autosavePeriod          = 100;
  
  //Modified Counter
  Bool_t   f_nparticle               = kFALSE;
  Bool_t   f_randomseed              = kFALSE;
  Bool_t   f_outputFilename          = kFALSE;
  Bool_t   f_startingMinimumDistance = kFALSE;
  Bool_t   f_derivingDegreeCriterion = kFALSE;
  Bool_t   f_derivingMinimumPoint    = kFALSE;
  Bool_t   f_velocity                = kFALSE;
  Bool_t   f_vperiod                 = kFALSE;
  Bool_t   f_autosavePeriod          = kFALSE;

  while( (c = getopt(argc, argv, "n:r:o:d:D:p:V:v:S:h")) != -1){
  // -1 means getopt() parse all options
  switch(c){
    case 'n':
      nparticle = std::stoi(optarg);
      f_nparticle = kTRUE;
      break;
    case 'r':
      randomseed = std::stoi(optarg);
      f_randomseed = kTRUE;
      break;
    case 'o':
      outputFilename = TString(optarg);
      f_outputFilename = kTRUE;
      break;
    case 'd':
      derivingDegreeCriterion = std::stod(optarg);
      f_derivingDegreeCriterion = kTRUE;
      break;
    case 'D':
      startingMinimumDistance = std::stod(optarg);
      f_startingMinimumDistance = kTRUE;
      break;
    case 'p':
      derivingMinimumPoint = std::stoi(optarg);
      f_derivingMinimumPoint = kTRUE;
      break;
    case 'V':
      velocity = std::stod(optarg);
      f_velocity = kTRUE;
      break;
    case 'v':
      vperiod = std::stod(optarg);
      f_vperiod = kTRUE;
      break;
    case 'S':
      autosavePeriod = std::stoi(optarg);
      f_autosavePeriod = kTRUE;
      break;
    case 'h':
      std::cout<<"Rutherford Experiment Simulation (Numerical)"<<std::endl;
      std::cout<<"If there's less argv than parameters, it will use default value for absents"<<std::endl<<std::endl;
      std::cout<<"-n \tNumber of MonteCarlo Simulation" <<std::endl;
      std::cout<<"-r \tRandom Seed (UInt)" <<std::endl;
      std::cout<<"-o \tOutput Filename (path)" <<std::endl;
      std::cout<<"-d \tAsymptote Degree Criterion" <<std::endl;
      std::cout<<"-D \tMinimum Incident Distanc" <<std::endl;
      std::cout<<"-p \tMinumum Number of Derived Points" <<std::endl;
      std::cout<<"-V \tIncident Particle's Velocity (r.SOL)" <<std::endl;
      std::cout<<"-v \tVerbose Period" <<std::endl;
      std::cout<<"-S \tAutosave Period" <<std::endl;
      std::cout<<"-h \tHelp (this message)" <<std::endl;
      return 0;
    case '?':
      printf("Unknown flag : %c", optopt);
      break;
    }
  }

  std::cout<<"===========================================" << std::endl;
  std::cout<<"Rutherford Experiment Simulation" << std::endl;
  std::cout<<"with following parameters" << std::endl << std::endl;

  std::cout<<"Number of Particles (n) : "<<  nparticle <<" " << (!f_nparticle ? "(Default)" : " ") << std::endl;
  std::cout<<"Random Seed (r) : "<< randomseed <<" " << (!f_randomseed ? "(Default)" : " ") << std::endl;

  std::cout<<"Asymptote Degree Criterion (d) : " << derivingDegreeCriterion << " " << (!f_derivingDegreeCriterion ? "(Default)" : " ") << std::endl;
  std::cout<<"Minumum Number of Derived Points (p) : " << derivingMinimumPoint << " " << (!f_derivingMinimumPoint ? "(Default)" : " ") << std::endl;
  std::cout<<"Incident Particle's Velocity (r.SOL) (V) : " << velocity << " " << (!f_velocity ? "(Default)" : " ") << std::endl;
  std::cout<<"Minimum Incident Distance (D) : " << startingMinimumDistance << " " << (!f_startingMinimumDistance ? "(Default)" : " ") << std::endl;

  std::cout<<"Output File's Name (o) : "<< outputFilename << " " << (!f_nparticle ? "(Default)" : " ") << std::endl;
  std::cout<<"Verbose Period (v) : " << vperiod << " " << (!f_vperiod ? "(Default)" : " ") << std::endl;
  std::cout<<"Autosave Period (S) : " << autosavePeriod << " " << (!f_autosavePeriod ? "(Default)" : " ") << std::endl;
  std::cout<<"===========================================" << std::endl << std::endl;

  return Experiment(nparticle, randomseed, outputFilename, derivingDegreeCriterion, derivingMinimumPoint, velocity, startingMinimumDistance, vperiod, autosavePeriod);

}
