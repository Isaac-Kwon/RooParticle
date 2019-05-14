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

int main(){
  event * eventT = new event();

  coulombForce * cp = new coulombForce();
  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  Double_t x1_[] = {0., 0., 0.};
  Double_t v1_[] = {0., 0., 0.};
  TVectorD x1 = TVectorD(3, x1_);
  TVectorD v1 = TVectorD(3, v1_);
  EMparticle * p1 = new EMparticle(206,79, x1, v1, true, false);
  eventT->AddParticle(p1);

  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  FixedTargetExperiment * FTE = new FixedTargetExperiment();
  FTE->setTemplate(eventT);

  TVectorD x2 = TVectorD(3);
  TVectorD v2 = TVectorD(3);

  Float_t imp_min = 0.;
  Float_t imp_max = +1000.;
  Float_t imp;

  EMparticle * p2;

  Double_t v2_[] = {0.05,0.,0.};
  v2 = TVectorD(3,v2_);

  Double_t x2_[] = {-200., 0., 0.,};

  Int_t i;
  Int_t j=0;

  for(j=0;j<1;j++){
    TFile * hfile = new TFile(TString::Format("Data/test_%d.root",j),"RECREATE");
    TTree * tree = new TTree(TString::Format("TT_%d",j), "Rutherford Scattering Angle");
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
    tree->Branch("SAngle", & SAngle_tree, "angle/D");



    for(i=0; i<1 ;i++ ){
    for(i=0; i<1000000 ;i++ ){
      imp = RandomFloat(imp_min,imp_max);
      x2_[1] = imp;
      x2 = TVectorD(3,x2_);
      p2 = new EMparticle(4,2, x2, v2, false, true);
      FTE->makeEvent(p2);
      FTE->getEvent()->DeriveMAX();
      imp_tree = imp;
      NPOINT_tree = p2->GetPath()->GetMaxNumber();
      vxF_tree = p2->GetPath()->GetLastV().operator[](0);
      vyF_tree = p2->GetPath()->GetLastV().operator[](1);
      SAngle_tree = angleXD(p2->GetPath()->GetLastV());
      tree->Fill();

      if(i%1000==0){
        cout<<j<<"/"<<i<<":"<<imp_tree<<":"<<vxF_tree<<":"<<SAngle_tree<<endl;
      }

      cout<<p1->GetQ()<<endl;
      cout<<p2->GetQ()<<endl;
      cout<<p1->GetM()<<endl;
      cout<<p2->GetM()<<endl;
      FTE->offEvent();
      FTE->delEvent();
    }
    tree->Print();
    tree->AutoSave();
    hfile->Write();
    hfile->Close();
  }






  return 0;

}
