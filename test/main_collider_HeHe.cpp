#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"

#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

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

int main(){

  event * eventT = new event();

  coulombForce * cp = new coulombForce();
  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  Double_t x1_[] = {1000., 0., 0.};
  Double_t v1_[] = {-0.05,0.,0.};
  TVectorD x1 = TVectorD(3, x1_);
  TVectorD v1 = TVectorD(3, v1_);

  // eventT->AddParticle(p1);
  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  ColliderExperiment * CDE = new ColliderExperiment();
  CDE->setTemplate(eventT);

  TVectorD x2 = TVectorD(3);
  TVectorD v2 = TVectorD(3);

  Float_t imp_min = 0.;
  Float_t imp_max = +1000.;
  Float_t imp;

  EMparticle * p1, * p2;

  Double_t v2_[] = {0.05,0.,0.};
  v2 = TVectorD(3,v2_);

  Double_t x2_[] = {-1000, 0., 0.,};

  Int_t i;
  Int_t j=0;

  for(j=0;j<1;j++){
    TFile * hfile = new TFile(TString::Format("Data/test_collider_long_%d.root",j),"RECREATE");
    TTree * tree = new TTree(TString::Format("T_%d",j), "Two Particle Colliding Angle");

    Double_t   imp_tree;
    Int_t      NPOINT_tree;
    Double_t   vx1_tree;
    Double_t   vx2_tree;
    Double_t   SAngle1_tree;
    Double_t   SAngle2_tree;
    Double_t   DCA_tree;
    Double_t   outtime_free;

    tree->Branch("imp", &imp_tree, "imp/D");
    tree->Branch("NPOINT", &NPOINT_tree, "Npoint/I");
    tree->Branch("vx1", & vx1_tree, "vx1/D");
    tree->Branch("vx2", & vx2_tree, "vx2/D");
    tree->Branch("SAngle1", & SAngle1_tree, "angle1/D");
    tree->Branch("SAngle2", & SAngle2_tree, "angle2/D");
    tree->Branch("DCA", & DCA_tree, "DCA/D");



    for(i=0; i<10000 ;i++ ){
      imp = RandomFloat(imp_min,imp_max);
      x2_[1] = imp;
      x2 = TVectorD(3,x2_);

      p1 = new EMparticle(4,2, x1, v1, false, true);
      p2 = new EMparticle(4,2, x2, v2, false, true);

      // p1->GetX().Print();
      // p2->GetX().Print();

      CDE->makeEvent(p2);
      CDE->getEvent()->AddParticle(p1);

      CDE->getEvent()->DeriveMAX(1);

      imp_tree = imp;
      NPOINT_tree = p2->GetPath()->GetMaxNumber();
      SAngle1_tree = angleXD(p2->GetPath()->GetLastV());
      SAngle2_tree = angleXD(CDE->getEvent()->getParticle(1)->GetPath()->GetLastV());
      p2->GetPath()->GetDCA(p1->GetPath(), DCA_tree, outtime_free);
      vx1_tree = p1->GetPath()->GetLastX()[0];
      vx2_tree = p2->GetPath()->GetLastX()[0];

      tree->Fill();

      if(i%100==0){
        cout<<j<<"/"<<i<<":"<<imp_tree<<":"<<SAngle1_tree<<":"<<SAngle2_tree<<":"<<DCA_tree<<":"<<outtime_free<<endl;
        cout<<"LAST POS P1: "<<p1->GetPath()->GetLastX()[0]<<endl;
        cout<<"LAST POS P2: "<<p2->GetPath()->GetLastX()[0]<<endl;
        // p2->GetPath()->GetLastV().Print();
        // CDE->getEvent()->getParticle(1)->GetPath()->GetLastV().Print();
      }
      CDE->offEvent();
      CDE->delEvent();

    }
    tree->Print();
    tree->AutoSave();
    // hfile->Write();
    hfile->Close();
    delete hfile;
  }

  return 0;
}
