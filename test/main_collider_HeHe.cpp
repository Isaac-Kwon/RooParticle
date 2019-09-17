#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"

#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

class ColliderExperiment{
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

void ColliderExperiment::makeEvent(particle* p){
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
  return result*(180./3.1415926535);
}

int main(int argc, char** argv){

  event * eventT = new event();

  coulombForce * cp = new coulombForce();
  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  Double_t x1_[] = {1000., 0., 0.};
  Double_t v1_[] = {-0.05,0.,0.};
  TVector3 x1 = TVector3(x1_);
  TVector3 v1 = TVector3(v1_);

  // eventT->AddParticle(p1);
  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  ColliderExperiment * CDE = new ColliderExperiment();
  CDE->setTemplate(eventT);

  TVector3 x2 = TVector3();
  TVector3 v2 = TVector3();

  Float_t imp_min = 0.;
  Float_t imp_max = +1000.;

  EMparticle * p1, * p2;

  Double_t v2_[] = {0.05,0.,0.};
  v2 = TVector3(v2_);

  Double_t x2_[] = {-1000, 0., 0.,};

  UInt_t randseed = (unsigned int) time(NULL);
  TRandom3 * rg = new TRandom3(randseed);

  Int_t i;
  Int_t j=0;

  TFile * hfile;

  
  TString findex;
  if(argc==0){
    findex = TString("0");
  }else{
    findex = TString(argv[1]);
  }
  
  hfile = new TFile(TString::Format("Data/test_collider_%s.root",argv[1]),"RECREATE");
  TTree * tree = new TTree("T_Colliding", "Two Particle Colliding Angle");

  Float_t    imp;
  Int_t      NPOINT;
  Double_t   vx1;
  Double_t   vx2;
  Double_t   SAngle1;
  Double_t   SAngle2;
  Double_t   DCA;
  Double_t   outtime;

  tree->Branch("imp", &imp, "imp/F");
  tree->Branch("NPOINT", &NPOINT, "Npoint/I");
  tree->Branch("x1", & vx1, "x1/D");
  tree->Branch("x2", & vx2, "x2/D");
  // tree->Branch("x1p", vx1p, "x1p[Npoint]/D");
  // tree->Branch("x2p", vx1p, "x2p[Npoint]/D");
  tree->Branch("SAngle1", & SAngle1, "angle1/D");
  tree->Branch("SAngle2", & SAngle2, "angle2/D");
  tree->Branch("DCA", & DCA, "DCA/D");

  cout<<"Starting Data Production Index "<<findex<<endl;
  cout<<"RandomSeed : "<<randseed<<endl;
  for(i=0; i<10000; i++){
    imp = rg->Uniform(imp_min,imp_max);
    x2_[1] = imp;
    x2 = TVector3(x2_);

    p1 = new EMparticle(4,2, x1, v1, false, true);
    p2 = new EMparticle(4,2, x2, v2, false, true);

    CDE->makeEvent(p2);
    CDE->getEvent()->AddParticle(p1);
    CDE->getEvent()->DeriveMAX(1);

    NPOINT = p2->GetPath()->GetMaxNumber();
    SAngle1 = angleXD(p2->GetPath()->GetLastV());
    SAngle2 = angleXD(CDE->getEvent()->getParticle(1)->GetPath()->GetLastV());
    p2->GetPath()->GetDCA(p1->GetPath(), DCA, outtime);
    vx1 = p1->GetPath()->GetLastV()[0];
    vx2 = p2->GetPath()->GetLastV()[0];

    tree->Fill();

    if(i%100==0){
      cout<<findex<<":"<<"/"<<i<<":"<<imp<<":"<<SAngle1<<":"<<SAngle2<<":"<<DCA<<":"<<outtime<<endl;
      cout<<"LAST POS P1: "<<p1->GetPath()->GetLastX()[0]<<endl;
      cout<<"LAST POS P2: "<<p2->GetPath()->GetLastX()[0]<<endl;
      tree->AutoSave();
    }
      // p2->GetPath()->GetLastV().Print();
      // CDE->getEvent()->getParticle(1)->GetPath()->GetLastV().Print();
    CDE->offEvent();
    CDE->delEvent();
  }


  if(argc==0){
    tree->Print();
  }else{
    cout<<"JOB END: "<<argv[1]<<endl;
  }
  tree->AutoSave();
  hfile->Write();
  hfile->Close();
  delete hfile;


  return 0;
}
