#include "TROOT.h"
#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"

#include "iostream"
#include <math.h>
#include <ctime>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"

//==================================================================//
// 3 body experiment with non-uniform process to source plain
// Monte-carlo parameter
// "2 of Impact Parameters -bA, bB- and Projectile Angle to Target -Xi-"
// for all parameter, TRandom3 with seed of initial simulating time is used.
//
// Fixed Variable
// 1. Target Position (0,0,0) [fm] : fixed eternally
// 2. Projectile speed (0.05) [fm/c] : near Ek = 4.66MeV for alpha
// 3. ProjectileA azimuthal angle (theta = 0)
// 
// Independent Variable :: Randomly generated with TRandom3
// 1. bA: Impact Parameter for Projectile A to Fixed Target [0,1000] [fm] 
// 2. bB: Impact Parameter for Projectile B to Fixed Target [0,1000] [fm]
// 3. xi: Angle of Projected Impact Paramter Angle of ProjectileA to Projected Impact Paramter Angle of ProjectileB [0,2*pi] [rad]
//
// Dependent Variable :: Calculated from each motion-characteristics
//
// Initally Selected Event (filtered right after random-number generation
// Nan
//  
//==================================================================//



class MBodyExperiment{
public:
  void setTemplate(event* eventT_){eventT=eventT_;}
  void makeEvent(particle* p1, particle* p2);
  void offEvent(){eventOn=false;}
  void delEvent(){eventOn=false; delete currentEvent;}
  event* getEvent(){return currentEvent;}
protected:
private:
  event * eventT;
  Bool_t eventOn = false;
  event * currentEvent;
};

void MBodyExperiment::makeEvent(particle* p1, particle* p2){
  eventOn = true;
  currentEvent = new event(*eventT);
  currentEvent->AddParticle(p1);
  currentEvent->AddParticle(p2);
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

// Multifile inplemented
// Arguments
// 1. file suffix (default: 0)
// Running Example 1: ./test_3body test_1
// Running Example 2: ./test_3body test_02
// Running Example 2: ./test_3body 2

int main(int argc, char** argv){

  //event template definition
  event * eventT = new event();

  //base volume and force
  coulombForce * cp = new coulombForce();
  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  //fixed target definition
  Double_t x1_[] = {0., 0., 0.};
  Double_t v1_[] = {0., 0., 0.};
  TVectorD x1 = TVectorD(3, x1_);
  TVectorD v1 = TVectorD(3, v1_);
  EMparticle * p1 = new EMparticle(196,79, x1, v1, true, false);

  //event template
  eventT->AddParticle(p1);
  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  //**************************//
  //event template end
  //**************************//

  //Multi body experiment definition
  MBodyExperiment * MBE = new MBodyExperiment();
  MBE->setTemplate(eventT); //template input

  TVectorD x2 = TVectorD(3);
  TVectorD v2 = TVectorD(3);
  TVectorD x3 = TVectorD(3);
  TVectorD v3 = TVectorD(3);

  //Select Criterion to Random Variables.
  Float_t ImpactAMin = 0.;
  Float_t ImpactAMax = 1000.;
  Float_t ImpactBMin = 0.;
  Float_t ImpactBMax = 1000.;

  Float_t ImpactAngleMin = 0.;
  Float_t ImpactAngleMax = TMath::Pi();

  //Starting Z
  Double_t startAZ =  1000;
  Double_t startBZ = -1000;

  //particle Definition

  EMparticle *p2, *p3;

  Double_t x2_[] = { 0., 0., startAZ};
  Double_t x3_[] = { 0., 0., startBZ};

  Double_t v2_[] = { 0., 0. , -0.05 };
  Double_t v3_[] = { 0., 0. , 0.05 };
  v2 = TVectorD(3,v2_);
  v3 = TVectorD(3,v3_);

  //Tree data construction
  TString rfilename;
  if(argc==1){
    rfilename = TString("Data/MultiBody_Modified_0.root");
  }else if(argc==2){
    rfilename = TString::Format("Data/MultiBody_Modified_%s.root", argv[1]);
  }

  UInt_t randseed = (unsigned int) time(NULL);

  cout<<"SAVE DATA IN: "<<rfilename<<std::endl;
  cout<<"Random Seed: "<<randseed<<std::endl;
  TFile * hfile = new TFile(rfilename,"RECREATE");
  TTree * tree = new TTree("T_Multibody", "Multibody particle interaction");


  //starting point 1
  Double_t   sx1;
  Double_t   sy1;
  Double_t   sz1;
  Double_t   ImpactA; //Impact Paramter A to Fixed Target 

  //starting point 2
  Double_t   sx2;
  Double_t   sy2;
  Double_t   sz2;
  Double_t   ImpactB; //Impact Paramter B to Fixed Target 

  //Impact Angle
  Double_t   ImpactAngle;

  //last velocity 1
  Double_t   lvx1;
  Double_t   lvy1;
  Double_t   lvz1;

  //last velocity 2
  Double_t   lvx2;
  Double_t   lvy2;
  Double_t   lvz2;

  //Distance Closest Approach
  Double_t   DCATA;
  Double_t   DCATB;
  Double_t   DCAAB;

  //Time Closest Approach
  Double_t   TCATA;
  Double_t   TCATB;
  Double_t   TCAAB;

  // Double_t   outtime;
  Int_t      RN;

  //Branch: starting point 1
  tree->Branch("sx1", &sx1 , "sx1/D");
  tree->Branch("sy1", &sy1 , "sy1/D");
  tree->Branch("sz1", &sz1 , "sz1/D");
  tree->Branch("Impact1", &ImpactA , "Impact1/D");

  //Branch: starting point 2
  tree->Branch("sx2", &sx2 , "sx2/D");
  tree->Branch("sy2", &sy2 , "sy2/D");
  tree->Branch("sz2", &sz2 , "sz2/D");
  tree->Branch("Impact2", &ImpactB , "Impact2/D");

  //Branch: ImpactAngle
  tree->Branch("ImpactAngle", &ImpactAngle , "ImpactAngle/D");

  //Branch: last velocity 1
  tree->Branch("lvx1", &lvx1, "lvx1/D");
  tree->Branch("lvy1", &lvy1, "lvy1/D");
  tree->Branch("lvz1", &lvz1, "lvz1/D");

  //Branch: last velocity 2
  tree->Branch("lvx2", &lvx2, "lvx2/D");
  tree->Branch("lvy2", &lvy2, "lvy2/D");
  tree->Branch("lvz2", &lvz2, "lvz2/D");

  //Branch: DCA (TA, TB, AB)
  tree->Branch("DCA_TA", &DCATA, "DCA_TA/D");
  tree->Branch("DCA_TB", &DCATB, "DCA_TB/D");
  tree->Branch("DCA_AB", &DCAAB, "DCA_AB/D");

  //Branch: TCA (TA, TB, AB)
  tree->Branch("TCA_TA", &TCATA, "TCA_TA/D");
  tree->Branch("TCA_TB", &TCATB, "TCA_TB/D");
  tree->Branch("TCA_AB", &TCAAB, "TCA_AB/D");

  tree->Branch("RSeed",  &randseed,  "RSeed/i");
  tree->Branch("RN",  &RN,  "RN/I");

  TRandom3 * rg = new TRandom3(randseed);

  Int_t i;
  for(i=0; i<10000 ;i++ ){
    RN = i;
    ImpactA = rg->Uniform(ImpactAMin, ImpactAMax);
    ImpactB = rg->Uniform(ImpactBMin, ImpactBMax);
    ImpactAngle = rg->Uniform(ImpactAngleMin, ImpactAngleMax);

    x2_[0] = ImpactA;
    x2 = TVectorD(3,x2_);
    p2 = new EMparticle(4,2, x2, v2, false, true);

    x3_[0] = ImpactB*TMath::Cos(ImpactAngle);
    x3_[1] = ImpactB*TMath::Sin(ImpactAngle);
    x3 = TVectorD(3,x3_);
    p3 = new EMparticle(4,2, x3, v3, false, true);

    MBE->makeEvent(p2, p3);
    MBE->getEvent()->DeriveMAX();

    sx1  = x2_[0];
    sy1  = x2_[1];
    sz1  = x2_[2];

    sx2  = x3_[0];
    sy2  = x3_[1];
    sz2  = x3_[2];

    lvx1 = p2->GetPath()->GetLastV().operator[](0);
    lvy1 = p2->GetPath()->GetLastV().operator[](1);
    lvz1 = p2->GetPath()->GetLastV().operator[](2);

    lvx2 = p3->GetPath()->GetLastV().operator[](0);
    lvy2 = p3->GetPath()->GetLastV().operator[](1);
    lvz2 = p3->GetPath()->GetLastV().operator[](2);

    p2->GetPath()->GetDCA(p3->GetPath(), DCAAB, TCAAB);
    p2->GetPath()->GetDCA(x1, DCATA, TCATA);
    p3->GetPath()->GetDCA(x1, DCATB, TCATB);

    tree->Fill();

    if(i%100==0){
      cout<<(argc>2?argv[2]:"0")<<"/"<<i<<" : DCA "<<DCAAB<<"|"<<DCATA<<","<<DCATB<<" : TCA "<<TCAAB<<"|"<<TCATA<<","<<TCATB<<endl;
    }
    tree->AutoSave();

    MBE->offEvent();
    MBE->delEvent();
  }

  tree->Print();
  cout<<"END PRODUCTION "<<(argc>2?argv[2]:"")<<endl;
  tree->AutoSave();
  hfile->Close();
  delete hfile;
  // tree->Delete();
  // cout<<"TREE DELETED"<<endl;

  return 0;
}
