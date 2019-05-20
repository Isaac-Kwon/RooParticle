#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"

#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"


using namespace std;

// Test::Rutherford Scattering Tress Analysis

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


Float_t RandomFloat(Float_t a, Float_t b) {
    Float_t random = ((Float_t) rand()) / (Float_t) RAND_MAX;
    Float_t diff = b - a;
    Float_t r = random * diff;
    return a + r;
}

TVectorD unitv(TVectorD v){
  TVectorD result = TVectorD(3);
  result = v*(1/sqrt(v.Norm2Sqr()));
  return result;
}

Double_t angleXD(TVectorD v){
  TVectorD unit = TVectorD(3);
  unit = unitv(v);
  Double_t result = acos(unit[0]);
  return result*(180./3.1415926535);
}


//Multifile inplemented
//Arguments
//1. exp name (default: test)
//2. file number (default: 0)
// Running Example 1: ./main 1         :: file will be FixedTarget79_2_4_1.root
// Running Example 2: ./main  02       :: file will be FixedTarget79_2_4_02.root

int main(int argc, char** argv){

  /* START OF EVENT Template Definition */
  /* Event Template */
  // TEMPLATE PARTICLE
  // 1. Fixed Target GOLD
  //   - Atomic Number: 79
  //   - Mass Number 197
  //   - At position (0,0,0)
  //   - with initial Velocity (0,0,0)
  //   - Rigid Body (not affected by any force. no momemtum variation)
  // TEMPLATE FORCE
  // 1. Coulomb Force (in force.cpp/hpp)
  // TEMPLATE VOLUME
  // (There is volume def. but don't do anything)

  event * eventT = new event();
  coulombForce * cp = new coulombForce();
  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  Double_t x1_[] = {0., 0., 0.};
  Double_t v1_[] = {0., 0., 0.};
  TVectorD x1 = TVectorD(3, x1_);
  TVectorD v1 = TVectorD(3, v1_);
  EMparticle * p1 = new EMparticle(197,79, x1, v1, true, false);

  eventT->AddParticle(p1);
  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  /* END OF EVENT Template Definition */
  /* Fixed Target Experiment Production */

  FixedTargetExperiment * FTE = new FixedTargetExperiment();

  FTE->setTemplate(eventT); // Set Event Tempate with defined template

  /* START OF PROJECTILE PARTICLE DEFINITION */
  TVectorD x2 = TVectorD(3); // DEF) Projectile Particle Initial Position (Vector)
  TVectorD v2 = TVectorD(3); // DEF) Projectile Particle Initial Velocity (Vector)

  Float_t imp_min = 0.;      // DEF,SET) Projectile Particle Minimum Impact Parameter
  Float_t imp_max = +1000.;  // DEF,SET) Projectile Particle Maximum Impact Parameter
  Float_t imp;               // DEF)     Projectile Particle Impact Parameter

  EMparticle * p2;           // DEF)     Projectile Particle

  Double_t v2_[] = {0.05,0.,0.}; // DEF,SET) Projectile Particle Initial Velocity (NumArray)
  v2 = TVectorD(3,v2_);          // SET) Projectile Particle Initial Velocity (Vector)

  Double_t x2_[] = {-1000., 0., 0.,}; // DEF,SET) Projectile Particle Initial Position (NumArray)

  Int_t i; // Produced Particle Number
  Int_t j=0; // Produced ROOT File Number



  //Tree data construction
  //MultiTreading
  TString rfilename;

  if(argc==1){
    rfilename = TString("Data/FixedTarget79_2_4_0.root");
  }else if(argc==2){
    rfilename = TString::Format("Data/FixedTarget79_2_4_%s.root", argv[1]);
  }

  cout<<"SAVE DATA IN: "<<rfilename<<endl;
  TFile * hfile = new TFile(rfilename,"RECREATE");
  TTree * tree = new TTree("T_FixedGold", "Rutherford Scattering Angle");

  // TREE Entry Properties
  // for Each event, followings will be recorded.
  // imp_tree: Impact Parameter
  // NPOINT_tree: Number of points iterated.
  // vxF_tree: Final VelocityX
  // vyF_tree: Final VelocityY
  // vyF_tree: Scattering Angle (calculated from vxF, vyF)

  Double_t   imp_tree;
  Int_t      NPOINT_tree;
  Double_t   vxF_tree;
  Double_t   vyF_tree;
  Double_t   SAngle_tree;

  tree->Branch("imp", &imp_tree, "imp/D");
  tree->Branch("NPOINT", &NPOINT_tree, "Npoint/I");
  // tree->Branch("x", x_tree, "x[Npoint]/D");
  // tree->Branch("y", y_tree, "y[Npoint]/D");
  // tree->Branch("vx", vx_tree, "vx[Npoint]/D");
  // tree->Branch("vy", vy_tree, "vy[Npoint]/D");
  tree->Branch("vxF", & vxF_tree, "vxF/D");
  tree->Branch("vyF", & vyF_tree, "vyF/D");
  tree->Branch("SAngle", & SAngle_tree, "SAngle/D");


  for(i=0; i<10000; i++ ){
    //Start Single Event Production
    //1. for Randomized Impact Parameter
    imp = RandomFloat(imp_min,imp_max);
    x2_[1] = imp;
    x2 = TVectorD(3,x2_);
    //Produce Particle with impact parameter
    p2 = new EMparticle(4,2, x2, v2, false, true);
    //make Event with projectile p2
    FTE->makeEvent(p2);
    FTE->getEvent()->DeriveMAX();
    //take value from event
    imp_tree = imp;
    NPOINT_tree = p2->GetPath()->GetMaxNumber();
    vxF_tree = p2->GetPath()->GetLastV().operator[](0);
    vyF_tree = p2->GetPath()->GetLastV().operator[](1);
    SAngle_tree = angleXD(p2->GetPath()->GetLastV());
    tree->Fill();

    if(i%1000==0){
      cout<<j<<"/"<<i<<":"<<imp_tree<<":"<<vxF_tree<<":"<<SAngle_tree<<endl;
    }

    FTE->offEvent();
    FTE->delEvent();
    tree->AutoSave();
  }
  tree->Print();
  hfile->Write();
  hfile->Close();

  return 0;

}
