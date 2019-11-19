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
#include "recorder.hpp"

#include "iostream"
#include <math.h>
#include "unistd.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"

TApplication * app;

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

int Experiment(Double_t ImpactParameter,
               TString outputFilename,
               Double_t derivingDegreeCriterion, Int_t derivingMinimumPoint,
               Double_t velocity, Double_t startingMinimumDistance,
               Int_t autosavePeriod){
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

  Float_t imp;

  EMparticle * p2;

  Double_t v2_[] = {velocity,0.,0.};
  v2 = TVector3(v2_);

  Double_t x2_[] = {-5000., 0., 0.,};

  TString fname = outputFilename;

  TFile * hfile = new TFile(fname,"RECREATE");
  TTree * tree = new TTree("Track", "Rutherford Scattering Track");

  std::cout<< "ROOT File save in " << fname << std::endl;

  Double_t   imp_tree;
  Double_t   x_tree;
  Double_t   y_tree;
  Double_t   vx_tree;
  Double_t   vy_tree;
  Double_t   t_tree;
  

  Int_t DeriveN;


  tree->Branch("imp", &imp_tree, "imp/D");
  tree->Branch("DeriveN", &DeriveN, "DeriveN/I");
  tree->Branch("x",  & x_tree, "sx/D");
  tree->Branch("y",  & y_tree, "sy/D");
  tree->Branch("vx", & vx_tree, "vxF/D");
  tree->Branch("vy", & vy_tree, "vyF/D");
  tree->Branch("t", & t_tree, "t/D");

  Double_t ImpWall = startingMinimumDistance * TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());

    imp = ImpactParameter;

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
    FTE->getEvent()->AddInspector(new inspector(FTE->getEvent()->getParticle(0), p2, "SPDG", 0.04999));
    
    imp_tree = imp;
    for(DeriveN=0; !(FTE->getEvent()->Inspect()); DeriveN++){
      if(DeriveN%autosavePeriod==0) tree->AutoSave();
      x_tree = FTE->getEvent()->getParticle(1)->GetX().operator[](0);
      y_tree = FTE->getEvent()->getParticle(1)->GetX().operator[](1);
      vx_tree = FTE->getEvent()->getParticle(1)->GetV().operator[](0);
      vy_tree = FTE->getEvent()->getParticle(1)->GetV().operator[](1);
      t_tree  = FTE->getEvent()->getParticle(1)->GetTime();
      tree->Fill();
      FTE->getEvent()->DeriveDT();
      
    }

    x_tree = FTE->getEvent()->getParticle(1)->GetX().operator[](0);
    y_tree = FTE->getEvent()->getParticle(1)->GetX().operator[](1);
    vx_tree = FTE->getEvent()->getParticle(1)->GetV().operator[](0);
    vy_tree = FTE->getEvent()->getParticle(1)->GetV().operator[](1);
    t_tree  = FTE->getEvent()->getParticle(1)->GetTime();
    tree->Fill();

    FTE->offEvent();
    FTE->delEvent();

  tree->Print();

  TCanvas * c1 = new TCanvas("c1", "c1", 600, 400);
  c1->cd();
  tree->Draw("x:y", "", "goff");
  TGraph * g1 = new TGraph(tree->GetEntries(), tree->GetV1(), tree->GetV2());
  g1->Draw();
  
  app->Run();
  hfile->Close();
  delete hfile;

  return 0;
}

//===============================================================================
// ./test_rutherford_inspect_oneparticle_tracksave
// Rutherford Experiment Simulation (Numerical)
// If there's less argv than parameters, it will use default value for absents
// -b   Impact Parameter
// -o   Output Filename (path)
// -d   Asymptote Degree Criterion
// -D   Minimum Incident Distance
// -p   Minumum Number of Derived Points
// -V   Incident Particle's Velocity (r.SOL)
// -S   Autosave Period
// -h   Help
//===============================================================================

int main(int argc, char *argv[]){
  char c; // option

  //Default setting
  Double_t ImpactParameter         = 50.;
  TString  outputFilename          = "Data/test_rutherford_inspect_oneparticle_showtrack.root";
  Double_t derivingDegreeCriterion = 1;
  Int_t    derivingMinimumPoint    = 50000;
  Double_t startingMinimumDistance = 5000;
  Double_t velocity                = 0.05;
  Int_t    autosavePeriod          = 100;
  
  //Modified Counter
  Bool_t   f_impactparameter         = kFALSE;
  Bool_t   f_outputFilename          = kFALSE;
  Bool_t   f_derivingDegreeCriterion = kFALSE;
  Bool_t   f_derivingMinimumPoint    = kFALSE;
  Bool_t   f_startingMinimumDistance = kFALSE;
  Bool_t   f_velocity                = kFALSE;
  Bool_t   f_autosavePeriod          = kFALSE;

  while( (c = getopt(argc, argv, "b:o:d:D:p:V:S:h")) != -1){
  // -1 means getopt() parse all options
  switch(c){
    case 'b':
      ImpactParameter = std::stod(optarg);
      f_impactparameter = kTRUE;
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
    case 'S':
      autosavePeriod = std::stoi(optarg);
      f_autosavePeriod = kTRUE;
      break;
    case 'h':

    // -b   Impact Parameter
    // -o   Output Filename (path)
    // -d   Asymptote Degree Criterion
    // -D   Minimum Incident Distance
    // -p   Minumum Number of Derived Points
    // -V   Incident Particle's Velocity (r.SOL)
    // -S   Autosave Period
    // -h   Help
      std::cout<<"Rutherford Experiment Simulation (Numerical)"<<std::endl;
      std::cout<<"If there's less argv than parameters, it will use default value for absents"<<std::endl<<std::endl;
      std::cout<<"-b \tImpact Parameter" <<std::endl;
      std::cout<<"-o \tOutput Filename (path)" <<std::endl;
      std::cout<<"-d \tAsymptote Degree Criterion" <<std::endl;
      std::cout<<"-D \tMinimum Incident Distance" <<std::endl;
      std::cout<<"-p \tMinumum Number of Derived Points" <<std::endl;
      std::cout<<"-V \tIncident Particle's Velocity (r.SOL)" <<std::endl;
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

  std::cout<<"Impact Parameter (n) : "<<  ImpactParameter <<" " << (!f_impactparameter ? "(Default)" : " ") << std::endl;

  std::cout<<"Asymptote Degree Criterion (d) : " << derivingDegreeCriterion << " " << (!f_derivingDegreeCriterion ? "(Default)" : " ") << std::endl;
  std::cout<<"Minumum Number of Derived Points (p) : " << derivingMinimumPoint << " " << (!f_derivingMinimumPoint ? "(Default)" : " ") << std::endl;
  std::cout<<"Incident Particle's Velocity (r.SOL) (V) : " << velocity << " " << (!f_velocity ? "(Default)" : " ") << std::endl;
  std::cout<<"Minimum Incident Distance (D) : " << startingMinimumDistance << " " << (!f_startingMinimumDistance ? "(Default)" : " ") << std::endl;

  std::cout<<"Output File's Name (o) : "<< outputFilename << " " << (!f_outputFilename ? "(Default)" : " ") << std::endl;
  std::cout<<"Autosave Period (S) : " << autosavePeriod << " " << (!f_autosavePeriod ? "(Default)" : " ") << std::endl;
  std::cout<<"===========================================" << std::endl << std::endl;

  app = new TApplication("app",&argc, argv);
  // app->Run();

  return Experiment(ImpactParameter, outputFilename, derivingDegreeCriterion, derivingMinimumPoint, velocity, startingMinimumDistance, autosavePeriod);

}
