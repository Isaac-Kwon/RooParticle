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

class MultibodyExperiment{
public:
  void setTemplate(event* eventT_){eventT=eventT_;}
  void makeEvent(particle* p1, particle* p2); // if event already on,
  void offEvent(){eventOn=false;}
  void delEvent(){eventOn=false; delete currentEvent;}
  event* getEvent(){return currentEvent;}
protected:
private:
  event * eventT;
  Bool_t eventOn = false;
  event * currentEvent;
};

void MultibodyExperiment::makeEvent(particle* p1, particle* p2){
  eventOn = true;
  currentEvent = new event(*eventT);
  currentEvent->AddParticle(p1);
  currentEvent->AddParticle(p2);
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

int Experiment(Int_t nparticle, Int_t randomseed,
               TString outputFilename,
               Double_t derivingDegreeCriterion, Int_t derivingMinimumPoint,
               Double_t velocity, Double_t startingMinimumDistance,
               Double_t MinimumImpactParameter, Double_t MaximumImpactParameter,
               Int_t vperiod, Int_t autosavePeriod){
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

  MultibodyExperiment * Exp = new MultibodyExperiment();
  Exp->setTemplate(eventT);

  TVector3 x2 = TVector3();
  TVector3 v2 = TVector3();
  TVector3 x3 = TVector3();
  TVector3 v3 = TVector3();

  // Float_t imp_min = 0.;
  // Float_t imp_max = +3000.;
  Float_t imp;

  EMparticle *p2, *p3;

  Double_t v2_[] = {velocity,0.,0.};
  Double_t v3_[] = {-velocity,0.,0.};
  v2 = TVector3(v2_);
  v3 = TVector3(v3_);

  Double_t x2_[] = {-5000., 0., 0.,};
  Double_t x3_[] = {+5000., 0., 0.,};

  Int_t i;
  Int_t j=0;

  TString fname = outputFilename;

  TFile * hfile = new TFile(fname,"RECREATE");
  TTree * tree = new TTree(TString::Format("TT_%d",j), "Rutherford Scattering Angle");

  std::cout<< "ROOT File save in " << fname << std::endl;

  Int_t      NPOINT_tree;
  Double_t   fx1_tree, fx2_tree;
  Double_t   fy1_tree, fy2_tree;

  Double_t   sx1_tree, sx2_tree;
  Double_t   sy1_tree, sy2_tree;

  Double_t   vxF1_tree, vxF2_tree;
  Double_t   vyF1_tree, vyF2_tree;

  Double_t   SAngle1_tree, SAngle2_tree;
  Double_t   DCA01_tree, DCA02_tree, DCA12_tree;

  tree->Branch("NPOINT", & NPOINT_tree, "Npoint/I");
  tree->Branch("fx1", & fx1_tree, "fx1/D");
  tree->Branch("fy1", & fy1_tree, "fy1/D");
  tree->Branch("fx2", & fx2_tree, "fx2/D");
  tree->Branch("fy2", & fy2_tree, "fy2/D");

  tree->Branch("sx1", & sx1_tree, "sx1/D");
  tree->Branch("sy1", & sy1_tree, "sy1/D");
  tree->Branch("sx2", & sx2_tree, "sx2/D");
  tree->Branch("sy2", & sy2_tree, "sy2/D");

  tree->Branch("vxF1", & vxF1_tree, "vxF1/D");
  tree->Branch("vyF1", & vyF1_tree, "vyF1/D");
  tree->Branch("vxF2", & vxF2_tree, "vxF2/D");
  tree->Branch("vyF2", & vyF2_tree, "vyF2/D");

  tree->Branch("SAngle1", & SAngle1_tree, "Sangle1/D");
  tree->Branch("SAngle2", & SAngle2_tree, "Sangle2/D");

  tree->Branch("DCA01", & DCA01_tree, "DCA01/D");
  tree->Branch("DCA02", & DCA02_tree, "DCA02/D");
  tree->Branch("DCA12", & DCA12_tree, "DCA12/D");

  TRandom3 * r1 = new TRandom3(randomseed);

  Double_t ImpWall = startingMinimumDistance * TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());

  for(i=0; i<nparticle ;i++ ){

    imp = r1->Uniform(MinimumImpactParameter, MaximumImpactParameter);

    x2_[1] = imp;
    x3_[1] = imp;

    if(imp>ImpWall){
      x2_[0] = -1 * imp / TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());
      x3_[0] = +1 * imp / TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());
    }else{
      x2_[0] = -1 * startingMinimumDistance;
      x3_[0] = +1 * startingMinimumDistance;
    }

    x2 = TVector3(x2_);
    x3 = TVector3(x3_);
    p2 = new EMparticle(4,2, x2, v2, false, false);
    p3 = new EMparticle(4,2, x3, v3, false, false);

    Exp->makeEvent(p2, p3);

    Exp->getEvent()->AddInspector(new inspectorP(Exp->getEvent()->getParticle(0), p2, "DEG", derivingDegreeCriterion));
    Exp->getEvent()->AddInspector(new inspectorP(Exp->getEvent()->getParticle(0), p2, "CNT", derivingMinimumPoint));
    Exp->getEvent()->AddInspector(new inspectorP(Exp->getEvent()->getParticle(0), p3, "DEG", derivingDegreeCriterion));
    Exp->getEvent()->AddInspector(new inspectorP(Exp->getEvent()->getParticle(0), p3, "CNT", derivingMinimumPoint));

    Exp->getEvent()->AddInspector(new inspectorP(Exp->getEvent()->getParticle(0), p3, "SPDG", 0.049999));
    Exp->getEvent()->AddInspector(new inspectorP(Exp->getEvent()->getParticle(0), p3, "SPDG", 0.049999));

    
    recorderD * DCARecorder01 = new recorderD(numeric_limits<Double_t>::max());
    recorderD * DCARecorder02 = new recorderD(numeric_limits<Double_t>::max());
    recorderD * DCARecorder12 = new recorderD(numeric_limits<Double_t>::max());
    
    Int_t deriveN;
    Double_t distance01=numeric_limits<Double_t>::max();
    Double_t distance02=numeric_limits<Double_t>::max();
    Double_t distance12=numeric_limits<Double_t>::max();

    for(deriveN=0; !(Exp->getEvent()->Inspect()); deriveN++){
      Exp->getEvent()->DeriveDTN();
      distance01 = ((Exp->getEvent()->getParticle(0)->GetX())-(p2->GetX())).Mag();
      distance02 = ((Exp->getEvent()->getParticle(0)->GetX())-(p3->GetX())).Mag();
      distance12 = ((p3->GetX())-(p2->GetX())).Mag();
      DCARecorder01->MinSave(distance01);
      DCARecorder02->MinSave(distance02);
      DCARecorder12->MinSave(distance12);
    }

    NPOINT_tree = deriveN;
    
    sx1_tree = x2_[0];
    sy1_tree = x2_[1];

    sx2_tree = x3_[0];
    sy2_tree = x3_[1];
    
    fx1_tree = p2->GetX().operator[](0);
    fy1_tree = p2->GetX().operator[](1);
    fx2_tree = p3->GetX().operator[](0);
    fy2_tree = p3->GetX().operator[](1);

    vxF1_tree = p2->GetV().operator[](0);
    vyF1_tree = p2->GetV().operator[](1);
    vxF2_tree = p3->GetV().operator[](0);
    vyF2_tree = p3->GetV().operator[](1);

    SAngle1_tree = angleXD(p2->GetV());
    SAngle2_tree = angleXD(p3->GetV());

    DCA01_tree = DCARecorder01->GetData();
    DCA02_tree = DCARecorder02->GetData();
    DCA12_tree = DCARecorder12->GetData();

    tree->Fill();

    if(i%autosavePeriod==0){
      tree->AutoSave();
    }

    if(i%vperiod==0){
      std::cout<<j<<"/"<<i<<" | IMP: "<<imp<<" | xf_x: " << fx1_tree << " | vf_x : "<<vxF1_tree<<" | SAngle (DEG): "<<SAngle1_tree * TMath::RadToDeg() <<" | POINT: "<<NPOINT_tree<<endl;
      
    }
    Exp->offEvent();
    Exp->delEvent();

  }

  tree->Print();
  tree->AutoSave();
  hfile->Close();
  delete hfile;

  return 0;
}

//===============================================================================
// ./test_rutherford_recorder
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
  TString  outputFilename          = "Data/test_multibody_oneline.root";
  Double_t derivingDegreeCriterion = 1;
  Int_t    derivingMinimumPoint    = 50000;
  Double_t startingMinimumDistance = 5000;
  Double_t velocity                = 0.05;
  Int_t    vperiod                 = 5;
  Int_t    autosavePeriod          = 100;
  Double_t minIP                   = 0.;
  Double_t maxIP                   = 3000.;
  
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
  Bool_t   f_minIP                   = kFALSE;
  Bool_t   f_maxIP                   = kFALSE;

  while( (c = getopt(argc, argv, "n:r:o:d:D:p:V:v:S:m:M:h")) != -1){
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
    case 'm':
      minIP = std::stod(optarg);
      f_minIP = kTRUE;
      break;
    case 'M':
      maxIP = std::stod(optarg);
      f_maxIP = kTRUE;
      break;
    case 'h':
      std::cout<<"Rutherford Experiment Simulation (Numerical)"<<std::endl;
      std::cout<<"If there's less argv than parameters, it will use default value for absents"<<std::endl<<std::endl;
      std::cout<<"-n \tNumber of MonteCarlo Simulation" <<std::endl;
      std::cout<<"-r \tRandom Seed (UInt)" <<std::endl;
      std::cout<<"-o \tOutput Filename (path)" <<std::endl;
      std::cout<<"-d \tAsymptote Degree Criterion" <<std::endl;
      std::cout<<"-D \tMinimum Incident Distance" <<std::endl;
      std::cout<<"-p \tMinumum Number of Derived Points" <<std::endl;
      std::cout<<"-V \tIncident Particle's Velocity (r.SOL)" <<std::endl;
      std::cout<<"-m \tMinimum Impact parameter" <<std::endl;
      std::cout<<"-M \tMaximum Impact parameter" <<std::endl;
      // std::cout<<"-M \tMaximum Impact parameter" <<std::endl;
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

  std::cout<<"Minimum Impact Parameter (m): " << minIP << " " << (!f_minIP ? "(Default)" : " ") << std::endl;
  std::cout<<"Maximum Impact Parameter (m): " << maxIP << " " << (!f_maxIP ? "(Default)" : " ") << std::endl;

  std::cout<<"Output File's Name (o) : "<< outputFilename << " " << (!f_nparticle ? "(Default)" : " ") << std::endl;
  std::cout<<"Verbose Period (v) : " << vperiod << " " << (!f_vperiod ? "(Default)" : " ") << std::endl;
  std::cout<<"Autosave Period (S) : " << autosavePeriod << " " << (!f_autosavePeriod ? "(Default)" : " ") << std::endl;
  std::cout<<"===========================================" << std::endl << std::endl;

  return Experiment(nparticle, randomseed, outputFilename, derivingDegreeCriterion, derivingMinimumPoint, velocity, startingMinimumDistance, minIP, maxIP, vperiod, autosavePeriod);

}
