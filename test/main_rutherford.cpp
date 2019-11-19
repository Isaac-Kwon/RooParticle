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

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

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

int main(){
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
  Float_t imp_max = +1000.;
  Float_t imp;

  EMparticle * p2;

  Double_t v2_[] = {0.05,0.,0.};
  v2 = TVector3(v2_);

  Double_t x2_[] = {-1000., 0., 0.,};

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
    Double_t   DCA_tree;
    Double_t   outtime_free;

    tree->Branch("imp", &imp_tree, "imp/D");
    tree->Branch("NPOINT", &NPOINT_tree, "Npoint/I");
    tree->Branch("vxF", & vxF_tree, "vxF/D");
    tree->Branch("vyF", & vyF_tree, "vyF/D");
    tree->Branch("SAngle", & SAngle_tree, "angle/D");
    tree->Branch("DCA", & DCA_tree, "DCA/D");



    for(i=0; i<100 ;i++ ){
      imp = RandomFloat(imp_min,imp_max);
      x2_[1] = imp;
      x2 = TVector3(x2_);
      p2 = new EMparticle(4,2, x2, v2, false, true);

      FTE->makeEvent(p2);
      FTE->getEvent()->DeriveMAX();

      imp_tree = imp;
      NPOINT_tree = p2->GetPath()->GetMaxNumber();
      vxF_tree = p2->GetPath()->GetLastV().operator[](0);
      vyF_tree = p2->GetPath()->GetLastV().operator[](1);
      SAngle_tree = angleXD(p2->GetPath()->GetLastV());
      p2->GetPath()->GetDCA(x1, DCA_tree, outtime_free);

      tree->Fill();

      if(i%50==0){
        std::cout<<j<<"/"<<i<<":"<<imp_tree<<":"<<vxF_tree<<":"<<SAngle_tree<<":"<<DCA_tree<<":"<<outtime_free<<endl;
      }
      FTE->offEvent();
      FTE->delEvent();

    }
    tree->Print();
    tree->AutoSave();
    hfile->Write();
    hfile->Close();
    delete hfile;
  }
  return 0;
}
