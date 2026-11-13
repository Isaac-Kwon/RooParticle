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
#include "runinfo.hpp"

#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

class FixedTargetExperiment{
public:
  void setTemplate(oddevent* eventT_){eventT=eventT_;}
  void makeEvent(particle* p); // if event already on,
  void offEvent(){eventOn=false;}
  void delEvent(){eventOn=false; delete currentEvent;}
  oddevent* getEvent(){return currentEvent;}
protected:
private:
  oddevent * eventT;
  Bool_t eventOn = false;
  oddevent * currentEvent;
};

void FixedTargetExperiment::makeEvent(particle* p){
  eventOn = true;
  currentEvent = new oddevent(*eventT);
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
  return result;
}

// Usage: ./main [title] [jobIndex [seed]]
//   title     name of the run, selects the output file Data/<title>_<jobIndex>.root (default "test").
//             Taken as the title only if it is not a number, so that both
//             "./main 3" and the job runners (macro/runmulti.py : target title jobIndex) work.
//   jobIndex  index of this job (default 0)
//   seed      explicit random seed (default: runinfo::kDefaultSeed + jobIndex)
int main(int argc, char** argv){
  Int_t   a = 1;
  TString title = "test";
  if(argc > 1 && !TString(argv[1]).IsDigit()){
    title = argv[1];
    a = 2;
  }
  Int_t jobIndex = (argc > a) ? atoi(argv[a]) : 0;
  Long64_t seedArg = (argc > a+1) ? atoll(argv[a+1]) : -1;
  UInt_t randseed = runinfo::SeedFromJob(jobIndex, seedArg);

  oddevent * eventT = new oddevent();

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

  Double_t imp_min = 0.;
  Double_t imp_max = +1000.;
  Double_t imp;

  EMparticle * p2;

  Double_t v2_[] = {0.05,0.,0.};
  v2 = TVector3(v2_);

  Double_t x2_[] = {-1000., 0., 0.,};

  Int_t i;
  Int_t j=jobIndex;

  TString fname = TString::Format("Data/%s_%d.root",title.Data(),j);

  TFile * hfile = new TFile(fname,"RECREATE");
  TTree * tree = new TTree(TString::Format("TT_%d",j), "Rutherford Scattering Angle");

  std::cout<< "ROOT File save in " << fname << std::endl;

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

  TRandom3 * r1 = new TRandom3(randseed);
  std::cout<< "Random seed : " << randseed << " (job index " << jobIndex << ")" << std::endl;

  const Int_t nevents = 100;
  runinfo::WriteCommon(randseed, 1., eventT->GetIntegratorName().Data(), nevents,
                       "DEG<1 AND CNT>50000");

  for(i=0; i<nevents ;i++ ){

    imp = r1->Uniform(imp_min, imp_max);

    x2_[1] = imp;
    x2 = TVector3(x2_);
    p2 = new EMparticle(4,2, x2, v2, false, true);

    FTE->makeEvent(p2);

    FTE->getEvent()->AddInspector(new inspectorP(FTE->getEvent()->getParticle(0), p2, "DEG", 1));
    FTE->getEvent()->AddInspector(new inspectorP(FTE->getEvent()->getParticle(0), p2, "CNT", 50000));
    
    FTE->getEvent()->DeriveInspect(1);

    imp_tree = imp;
    NPOINT_tree = p2->GetPath()->GetMaxNumber();
    vxF_tree = p2->GetPath()->GetLastV().operator[](0);
    vyF_tree = p2->GetPath()->GetLastV().operator[](1);
    SAngle_tree = angleXD(p2->GetPath()->GetLastV());
    p2->GetPath()->GetDCA(x1, DCA_tree, outtime_free);

    tree->Fill();

    if(i%5==0){
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

  return 0;
}
