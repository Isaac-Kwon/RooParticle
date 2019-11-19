////////////////////////////////////////////
// 
// RooParticle Test Application
// (Fixed Target) Rutherford Experiment
// 
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "oddevent.hpp"
#include "particle.hpp"
#include "force.hpp"

#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

class CollidingExperiment{
public:
  void setTemplate(oddevent* eventT_){eventT=eventT_;}
  void makeEvent(particle* p1, particle* p2); // if event already on,
  void offEvent(){eventOn=false;}
  void delEvent(){eventOn=false; delete currentEvent;}
  oddevent* getEvent(){return currentEvent;}
protected:
private:
  oddevent * eventT;
  Bool_t eventOn = false;
  oddevent * currentEvent;
};

void CollidingExperiment::makeEvent(particle* p1, particle* p2){
  eventOn = true;
  currentEvent = new oddevent(*eventT);
  currentEvent->AddParticle(p1);
  currentEvent->AddParticle(p2);
}

Double_t angleXD(TVector3 v){
  return acos(v.Unit()[0]);
}

int main(){
  oddevent * eventT = new oddevent();

  coulombForce * cp = new coulombForce();
  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  CollidingExperiment * exp = new CollidingExperiment();
  exp->setTemplate(eventT);
  
  Float_t imp_min = 0.;
  Float_t imp_max = +100.;
  Float_t imp;

  EMparticle *p1, *p2;

  TVector3 x1, v1, x2, v2;
  Double_t x1_[] = {+1000., 0., 0.};
  Double_t v1_[] = {-0.05, 0., 0.};
  Double_t x2_[] = {-1000., 0., 0.,};
  Double_t v2_[] = {0.05,0.,0.};
  
  v1 = TVector3(v1_);
  v2 = TVector3(v2_);

  Int_t i;
  Int_t j=0;

  TString fname = TString::Format("Data/test_colliding_inspect_odd_%d.root",j);

  TFile * hfile = new TFile(fname,"RECREATE");
  TTree * tree = new TTree(TString::Format("TT_%d",j), "Rutherford Scattering Angle");

  std::cout<< "ROOT File save in " << fname << std::endl;

  Double_t   imp_tree;
  Int_t      NPOINT_tree;
  Double_t   vxF1_tree;
  Double_t   vyF1_tree;
  Double_t   SAngle1_tree;
  Double_t   vxF2_tree;
  Double_t   vyF2_tree;
  Double_t   SAngle2_tree;
  Double_t   DCA_tree;
  Double_t   outtime_free;


  tree->Branch("imp", &imp_tree, "imp/D");
  tree->Branch("NPOINT", &NPOINT_tree, "Npoint/I");
  tree->Branch("vxF1", & vxF1_tree, "vxF1/D");
  tree->Branch("vyF1", & vyF1_tree, "vyF1/D");
  tree->Branch("vx2F", & vxF2_tree, "vxF2/D");
  tree->Branch("vy2F", & vyF2_tree, "vyF2/D");
  tree->Branch("SAngle1", & SAngle1_tree, "angle1/D");
  tree->Branch("SAngle2", & SAngle2_tree, "angle2/D");
  tree->Branch("DCA", & DCA_tree, "DCA/D");

  TRandom3 * r1 = new TRandom3();

  for(i=0; i<1000 ;i++ ){

    imp = r1->Uniform(imp_min, imp_max);

    x1_[1] = imp/2;
    x2_[1] = -imp/2;
    x1 = TVector3(x1_);
    x2 = TVector3(x2_);

    p1 = new EMparticle(4,2, x1, v1, false, true);
    p2 = new EMparticle(4,2, x2, v2, false, true);

    exp->makeEvent(p1, p2);

    exp->getEvent()->AddInspector(new inspectorP(p1, p2, "DEG", 1));
    exp->getEvent()->AddInspector(new inspectorP(p1, p2, "CNT", 50000));
    
    exp->getEvent()->DeriveInspect(1);
    // FTE->getEvent()->DeriveDTN(1,10);

    imp_tree = imp;
    NPOINT_tree = p1->GetPath()->GetMaxNumber();
    // std::cout<<"1";
    vxF1_tree = p1->GetPath()->GetLastV().operator[](0);
    // std::cout<<"2";
    vyF1_tree = p1->GetPath()->GetLastV().operator[](1);
    // std::cout<<"3";
    SAngle1_tree = angleXD(p1->GetPath()->GetLastV());
    // std::cout<<"4";
    vxF2_tree = p2->GetPath()->GetLastV().operator[](0);
    // std::cout<<"5";
    vyF2_tree = p2->GetPath()->GetLastV().operator[](1);
    // std::cout<<"6";
    SAngle2_tree = angleXD(p2->GetPath()->GetLastV());
    // std::cout<<"7"<<std::endl;
    p2->GetPath()->GetDCA(x1, DCA_tree, outtime_free);

    tree->Fill();

    if(i%50==0){
      std::cout<<j<<"/"<<i<<":"<<imp_tree<<":"<<NPOINT_tree<<":"<<vxF1_tree<<":"<<SAngle1_tree<<":"<<DCA_tree<<":"<<outtime_free<<endl;
    }
    exp->offEvent();
    exp->delEvent();

  }

  // tree->Print();
  tree->AutoSave();
  hfile->Write();
  hfile->Close();
  delete hfile;

  return 0;
}
