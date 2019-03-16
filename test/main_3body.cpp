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

class MBodyExperiment{
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

//Multifile inplemented
//Arguments
//1. exp name (default: test)
//2. file number (default: 0)
// Running Example 1: ./test_3body test 1
// Running Example 2: ./test_3body test 02

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
  EMparticle * p1 = new EMparticle(206,79, x1, v1, true, false);

  //event template
  eventT->AddParticle(p1);
  eventT->AddForce(cp);
  eventT->SetVolume(vol);

  //
  //event template end
  //

  //Multi body experiment definition
  MBodyExperiment * MBE = new MBodyExperiment();
  MBE->setTemplate(eventT); //template input

  TVectorD x2 = TVectorD(3);
  TVectorD v2 = TVectorD(3);
  TVectorD x3 = TVectorD(3);
  TVectorD v3 = TVectorD(3);

  //Square source
  Float_t impactLX = 1000;
  Float_t impactLY = 1000;

  //z = -zmax plain
  Float_t xmin = -1 * impactLX;
  Float_t xmax = impactLX;
  Float_t ymin = -1 * impactLY;
  Float_t ymax = impactLY;

  //Starting Z
  Double_t startZl =  1000;
  Double_t startZr = -1000;

  //particle Definition

  EMparticle * p2, * p3;

  Double_t x2_[] = { 0., 0., startZl};
  Double_t x3_[] = { 0., 0., startZr};

  Double_t v2_[] = { 0., 0. , -0.05 };
  Double_t v3_[] = { 0., 0. , 0.05 };
  v2 = TVectorD(3,v2_);
  v3 = TVectorD(3,v3_);

  //Tree data construction

  TString rfilename;
  if(argc==1){
    rfilename = TString("Data/test_0.root");
  }else if(argc==2){
    rfilename = TString::Format("Data/%s_0.root", argv[1]);
  }else if(argc==3){
    rfilename = TString::Format("Data/%s_%s.root", argv[1], argv[2]);
  }else{
    rfilename = TString::Format("Data/%s_%s_%s.root", argv[1], argv[2], argv[3]);
  }

  cout<<"SAVE DATA IN: "<<rfilename<<endl;
  TFile * hfile = new TFile(rfilename,"RECREATE");
  TTree * tree = new TTree("T", "Multibody particle interaction");

  Int_t      NPOINT_tree;

  //starting point 1
  Double_t   sx1_tree;
  Double_t   sy1_tree;
  Double_t   sz1_tree;

  //starting point 2
  Double_t   sx2_tree;
  Double_t   sy2_tree;
  Double_t   sz2_tree;

  //last velocity 1
  Double_t   lvx1_tree;
  Double_t   lvy1_tree;
  Double_t   lvz1_tree;

  //last velocity 2
  Double_t   lvx2_tree;
  Double_t   lvy2_tree;
  Double_t   lvz2_tree;

  //Distance Closest Approach
  Double_t   DCATA_tree;
  Double_t   DCATB_tree;
  Double_t   DCAAB_tree;

  //Time Closest Approach
  Double_t   TCATA_tree;
  Double_t   TCATB_tree;
  Double_t   TCAAB_tree;

  // Double_t * x1p_tree;
  // Double_t * y1p_tree;
  // Double_t * z1p_tree;

  // Double_t * x2p_tree;
  // Double_t * y2p_tree;
  // Double_t * z2p_tree;

  Double_t   outtime_free;
  Int_t      RN_tree;

  unsigned int randseed = (unsigned int) time(NULL);

  if(argc>=4){
    randseed = randseed * atof(argv[3]) + atof(argv[2]);
  }

  //Branch: starting point 1
  tree->Branch("sx1", &sx1_tree , "sx1/D");
  tree->Branch("sy1", &sy1_tree , "sy1/D");
  tree->Branch("sz1", &sz1_tree , "sz1/D");

  //Branch: starting point 2
  tree->Branch("sx2", &sx2_tree , "sx2/D");
  tree->Branch("sy2", &sy2_tree , "sy2/D");
  tree->Branch("sz2", &sz2_tree , "sz2/D");

  //Branch: last velocity 1
  tree->Branch("lvx1", &lvx1_tree, "lvx1/D");
  tree->Branch("lvy1", &lvy1_tree, "lvy1/D");
  tree->Branch("lvz1", &lvz1_tree, "lvz1/D");

  //Branch: last velocity 2
  tree->Branch("lvx2", &lvx2_tree, "lvx2/D");
  tree->Branch("lvy2", &lvy2_tree, "lvy2/D");
  tree->Branch("lvz2", &lvz2_tree, "lvz2/D");

  //Branch: DCA (TA, TB, AB)
  tree->Branch("DCA_TA", &DCATA_tree, "DCA_TA/D");
  tree->Branch("DCA_TB", &DCATB_tree, "DCA_TB/D");
  tree->Branch("DCA_AB", &DCAAB_tree, "DCA_AB/D");

  //Branch: TCA (TA, TB, AB)
  tree->Branch("TCA_TA", &TCATA_tree, "TCA_TA/D");
  tree->Branch("TCA_TB", &TCATB_tree, "TCA_TB/D");
  tree->Branch("TCA_AB", &TCAAB_tree, "TCA_AB/D");

  tree->Branch("RSeed",  &randseed,  "RSeed/I");
  tree->Branch("RSeed",  &RN_tree,  "RN/I");

  TRandom3 * rg = new TRandom3(randseed);

  Int_t i;
  for(i=0; i<1000 ;i++ ){
    RN_tree = i;
    x2_[0] = rg->Uniform(xmin, xmax);
    x2_[1] = rg->Uniform(ymin, ymax);
    x2 = TVectorD(3,x2_);
    p2 = new EMparticle(4,2, x2, v2, false, true);

    x3_[0] = rg->Uniform(xmin, xmax);
    x3_[1] = rg->Uniform(ymin, ymax);
    x3 = TVectorD(3,x3_);
    p3 = new EMparticle(4,2, x3, v3, false, true);

    MBE->makeEvent(p2, p3);
    MBE->getEvent()->DeriveMAX();

    sx1_tree  = x2_[0];
    sy1_tree  = x2_[1];
    sz1_tree  = x2_[2];

    sx2_tree  = x3_[0];
    sy2_tree  = x3_[1];
    sz2_tree  = x3_[2];

    lvx1_tree = p2->GetPath()->GetLastV().operator[](0);
    lvy1_tree = p2->GetPath()->GetLastV().operator[](1);
    lvz1_tree = p2->GetPath()->GetLastV().operator[](2);

    lvx2_tree = p3->GetPath()->GetLastV().operator[](0);
    lvy2_tree = p3->GetPath()->GetLastV().operator[](1);
    lvz2_tree = p3->GetPath()->GetLastV().operator[](2);

    p2->GetPath()->GetDCA(p3->GetPath(), DCAAB_tree, TCAAB_tree);
    p2->GetPath()->GetDCA(x1, DCATA_tree, TCATA_tree);
    p3->GetPath()->GetDCA(x1, DCATB_tree, TCATB_tree);

    tree->Fill();

    if(i%100==0){
      cout<<(argc>2?argv[2]:"0")<<"/"<<i<<" : DCA "<<DCAAB_tree<<"|"<<DCATA_tree<<","<<DCATB_tree<<" : TCA "<<TCAAB_tree<<"|"<<TCATA_tree<<","<<TCATB_tree<<endl;
      tree->AutoSave();
    }
  }
  MBE->offEvent();
  MBE->delEvent();

  tree->Print();
  cout<<"END PRODUCTION "<<(argc>2?argv[2]:"")<<endl;
  tree->AutoSave();
  hfile->Close();
  delete hfile;
  // tree->Delete();
  // cout<<"TREE DELETED"<<endl;

  return 0;
}
