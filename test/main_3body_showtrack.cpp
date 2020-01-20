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
#include "experiment.hpp"

#include "iostream"
#include <math.h>
#include "unistd.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"

TApplication * app;

struct particleParameters{
  Double_t A;
  Double_t q;
  Double_t imp;
};

int Experiment(particleParameters pi1, particleParameters pi2,
               TString outputFilename,
               Double_t derivingDegreeCriterion, Int_t derivingMinimumPoint,
               Double_t velocity, Double_t startingMinimumDistance){

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

  experimentP * Exp = new experimentP(2);
  Exp->SetTemplate(eventT);

  experimentP * ExpF = new experimentP(1);
  ExpF->SetTemplate(eventT);

  TVector3 x2 = TVector3();
  TVector3 v2 = TVector3();
  TVector3 x3 = TVector3();
  TVector3 v3 = TVector3();

  EMparticle *p2, *p3, *pp;

  Double_t v2_[] = {velocity,0.,0.};
  Double_t v3_[] = {-velocity,0.,0.};
  v2 = TVector3(v2_);
  v3 = TVector3(v3_);

  Double_t x2_[] = {-5000., 0., 0.,};
  Double_t x3_[] = {+5000., 0., 0.,};

  TString fname = outputFilename;
  TFile * hfile = new TFile(fname,"RECREATE");
  TTree * tree = new TTree("TT", "Multibody Scattering Track");

  std::cout<< "ROOT File save in " << fname << std::endl;

  Int_t      NPOINT_tree;
  Double_t   t_tree;
  Double_t   x1_tree, x2_tree, xp_tree;
  Double_t   y1_tree, y2_tree, yp_tree;

  Double_t   vx1_tree, vx2_tree, vxp_tree;
  Double_t   vy1_tree, vy2_tree, vyp_tree;;

  tree->Branch("NPOINT", & NPOINT_tree, "Npoint/I");
  tree->Branch("t", &t_tree, "t/D");

  tree->Branch("x1", & x1_tree, "x1/D");
  tree->Branch("y1", & y1_tree, "y1/D");
  tree->Branch("x2", & x2_tree, "x2/D");
  tree->Branch("y2", & y2_tree, "y2/D");
  tree->Branch("xp", & xp_tree, "xp/D");
  tree->Branch("yp", & yp_tree, "yp/D");

  tree->Branch("vx1", & vx1_tree, "vx1/D");
  tree->Branch("vy1", & vy1_tree, "vy1/D");
  tree->Branch("vx2", & vx2_tree, "vx2/D");
  tree->Branch("vy2", & vy2_tree, "vy2/D");
  tree->Branch("vxp", & vxp_tree, "vxp/D");
  tree->Branch("vyp", & vyp_tree, "vyp/D");

  Double_t ImpWall = startingMinimumDistance * TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());

  x2_[1] = pi1.imp;
  x3_[1] = pi2.imp;

  Double_t impp = max(pi1.imp, pi2.imp);

  if(impp>ImpWall){
    x2_[0] = -1 * impp / TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());
    x3_[0] = +1 * impp / TMath::Tan(derivingDegreeCriterion*TMath::DegToRad());
  }else{
    x2_[0] = -1 * startingMinimumDistance;
    x3_[0] = +1 * startingMinimumDistance;
  }

  x2 = TVector3(x2_);
  x3 = TVector3(x3_);
  p2 = new EMparticle(pi1.A,pi1.q, x2, v2, false, false);
  p3 = new EMparticle(pi2.A,pi2.q, x3, v3, false, false);

  pp = new EMparticle(pi1.A,pi1.q, x2, v2, false, false);

  Exp->MakeEvent(p2, p3);
  ExpF->MakeEvent(pp);

  Exp->GetEvent()->AddInspector(new inspectorP(Exp->GetEvent()->getParticle(0), p2, "DEG", derivingDegreeCriterion));
  Exp->GetEvent()->AddInspector(new inspectorP(Exp->GetEvent()->getParticle(0), p2, "CNT", derivingMinimumPoint));
  Exp->GetEvent()->AddInspector(new inspectorP(Exp->GetEvent()->getParticle(0), p3, "DEG", derivingDegreeCriterion));
  Exp->GetEvent()->AddInspector(new inspectorP(Exp->GetEvent()->getParticle(0), p3, "CNT", derivingMinimumPoint));

  Exp->GetEvent()->AddInspector(new inspectorP(Exp->GetEvent()->getParticle(0), p3, "SPDG", 0.0499999));
  Exp->GetEvent()->AddInspector(new inspectorP(Exp->GetEvent()->getParticle(0), p3, "SPDG", 0.0499999));

  // ExpF->GetEvent()->AddInspector(new inspector(ExpF->GetEvent()->getParticle(0), pp, "SPDG", 0.050000));
  // ExpF->GetEvent()->AddInspector(new inspector(ExpF->GetEvent()->getParticle(0), pp, "DEG", derivingDegreeCriterion));
  // ExpF->GetEvent()->AddInspector(new inspector(ExpF->GetEvent()->getParticle(0), pp, "CNT", derivingMinimumPoint));
  
  Int_t deriveN;

  for(deriveN=0; !(Exp->GetEvent()->Inspect()); deriveN++){
    NPOINT_tree = deriveN;
    t_tree = Exp->GetEvent()->getParticle(1)->GetTime();
    x1_tree = Exp->GetEvent()->getParticle(1)->GetX().operator[](0);
    y1_tree = Exp->GetEvent()->getParticle(1)->GetX().operator[](1);
    x2_tree = Exp->GetEvent()->getParticle(2)->GetX().operator[](0);
    y2_tree = Exp->GetEvent()->getParticle(2)->GetX().operator[](1);
    xp_tree = ExpF->GetEvent()->getParticle(1)->GetX().operator[](0);
    yp_tree = ExpF->GetEvent()->getParticle(1)->GetX().operator[](1);

    vx1_tree = Exp->GetEvent()->getParticle(1)->GetV().operator[](0);
    vy1_tree = Exp->GetEvent()->getParticle(1)->GetV().operator[](1);
    vx2_tree = Exp->GetEvent()->getParticle(2)->GetV().operator[](0);
    vy2_tree = Exp->GetEvent()->getParticle(2)->GetV().operator[](1);
    vxp_tree = ExpF->GetEvent()->getParticle(1)->GetV().operator[](0);
    vyp_tree = ExpF->GetEvent()->getParticle(1)->GetV().operator[](1);
    tree->Fill();
    Exp->GetEvent()->DeriveDTN();
    ExpF->GetEvent()->DeriveDTN();
  }
  t_tree = Exp->GetEvent()->getParticle(1)->GetTime();
  x1_tree = Exp->GetEvent()->getParticle(1)->GetX().operator[](0);
  y1_tree = Exp->GetEvent()->getParticle(1)->GetX().operator[](1);
  x2_tree = Exp->GetEvent()->getParticle(2)->GetX().operator[](0);
  y2_tree = Exp->GetEvent()->getParticle(2)->GetX().operator[](1);
  vx1_tree = Exp->GetEvent()->getParticle(1)->GetV().operator[](0);
  vy1_tree = Exp->GetEvent()->getParticle(1)->GetV().operator[](1);
  vx2_tree = Exp->GetEvent()->getParticle(2)->GetV().operator[](0);
  vy2_tree = Exp->GetEvent()->getParticle(2)->GetV().operator[](1);
  vxp_tree = ExpF->GetEvent()->getParticle(1)->GetV().operator[](0);
  vyp_tree = ExpF->GetEvent()->getParticle(1)->GetV().operator[](1);
  tree->Fill();

  Exp->OffEvent();
  Exp->DelEvent();

  tree->Print();
  tree->AutoSave();

  TCanvas * c1 = new TCanvas("c1", "c1", 600, 400);
  TMultiGraph * mg = new TMultiGraph();

  c1->cd();

  tree->Draw("x1:y1","","goff");
  TGraph * g1 = new TGraph(tree->GetEntries(), tree->GetV1(), tree->GetV2());
  g1->SetMarkerColor(kRed);
  mg->Add(g1);

  tree->Draw("x2:y2","","goff");
  TGraph * g2 = new TGraph(tree->GetEntries(), tree->GetV1(), tree->GetV2());
  g2->SetMarkerColor(kBlue);
  mg->Add(g2);

  tree->Draw("xp:yp","","goff");
  TGraph * gp = new TGraph(tree->GetEntries(), tree->GetV1(), tree->GetV2());
  gp->SetMarkerColorAlpha(kBlack,0.8);
  mg->Add(gp);

  mg->Draw("AP");

  app->Run();
  hfile->Close();


  delete hfile;

  return 0;
}

int main(int argc, char *argv[]){
  char c; // option

  //Default setting
  Double_t q2                      = 2.;
  Double_t q3                      = 2.;
  Double_t A2                      = 4.;
  Double_t A3                      = 4.;
  Double_t imp1                    = 50;
  Double_t imp2                    = 50;
  TString  outputFilename          = "Data/test_multibody_showtrack.root";
  Double_t derivingDegreeCriterion = 1;
  Int_t    derivingMinimumPoint    = 50000;
  Double_t startingMinimumDistance = 5000;
  Double_t velocity                = 0.05;
  
  //Modified Counter
  Bool_t   f_AQ2                     = kFALSE;
  Bool_t   f_AQ3                     = kFALSE;
  Bool_t   f_imp1                    = kFALSE;
  Bool_t   f_imp2                    = kFALSE;
  Bool_t   f_outputFilename          = kFALSE;
  Bool_t   f_startingMinimumDistance = kFALSE;
  Bool_t   f_derivingDegreeCriterion = kFALSE;
  Bool_t   f_derivingMinimumPoint    = kFALSE;
  Bool_t   f_velocity                = kFALSE;

  while( (c = getopt(argc, argv, "q:Q:i:I:o:d:D:p:V:h")) != -1){
  // -1 means getopt() parse all options
  TString tempstr;
  TObjArray * temparr;
  switch(c){
    case 'q':
      // imp1 = std::stod(optarg);
      tempstr = TString(optarg);
      temparr = tempstr.Tokenize("/");
      A2=std::stod((((TObjString*) temparr->At(0))->String()).Data());
      q2=std::stod((((TObjString*) temparr->At(1))->String()).Data());
      f_AQ2 = kTRUE;
      break;
    case 'Q':
      tempstr = TString(optarg);
      temparr = tempstr.Tokenize("/");
      A3=std::stod((((TObjString*) temparr->At(0))->String()).Data());
      q3=std::stod((((TObjString*) temparr->At(1))->String()).Data());
      f_AQ3 = kTRUE;
      break;
    case 'i':
      imp1 = std::stod(optarg);
      f_imp1 = kTRUE;
      break;
    case 'I':
      imp2 = std::stod(optarg);
      f_imp2 = kTRUE;
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
    case 'h':
      std::cout<<"Rutherford Experiment Simulation (Numerical)"<<std::endl;
      std::cout<<"If there's less argv than parameters, it will use default value for absents"<<std::endl<<std::endl;
      std::cout<<"-q \t Form: A/Q, Mass Number and Charge (p2)" <<std::endl;
      std::cout<<"-Q \t Form: A/Q, Mass Number and Charge (p3)" <<std::endl;
      std::cout<<"-i \t(absolute) Impact Parameter (p2) (i)" <<std::endl;
      std::cout<<"-I \t(absolute) Impact Parameter (p3) (I)" <<std::endl;
      std::cout<<"-o \tOutput Filename (path)" <<std::endl;
      std::cout<<"-d \tAsymptote Degree Criterion" <<std::endl;
      std::cout<<"-D \tMinimum Incident Distance" <<std::endl;
      std::cout<<"-p \tMinumum Number of Derived Points" <<std::endl;
      std::cout<<"-V \tIncident Particle's Velocity (r.SOL)" <<std::endl;
      std::cout<<"-h \tHelp (this message)" <<std::endl;
      return 0;
    case '?':
      printf("Unknown flag : %c", optopt);
      break;
    }
  }

  particleParameters pi2 = particleParameters{A2, q2, imp1};
  particleParameters pi3 = particleParameters{A3, q3, imp2};

  std::cout<<"===========================================" << std::endl;
  std::cout<<"Rutherford Experiment Simulation" << std::endl;
  std::cout<<"with following parameters" << std::endl << std::endl;

  std::cout<<"(absolute) Impact Parameter (p2) (i) : "<< imp1 <<" " << (!f_imp1 ? "(Default)" : " ") << std::endl;
  std::cout<<"(absolute) Impact Parameter (p2) (i) : "<< imp1 <<" " << (!f_imp1 ? "(Default)" : " ") << std::endl;

  std::cout<<"Charge and Mass Number (p2) (q) : "<< q2 << "|" << A2 << " " << (!f_AQ2 ? "(Default)" : " ") << std::endl;
  std::cout<<"Charge and Mass Number (p3) (Q) : "<< q3 << "|" << A3 << " " << (!f_AQ3 ? "(Default)" : " ") << std::endl;

  std::cout<<"Asymptote Degree Criterion (d) : " << derivingDegreeCriterion << " " << (!f_derivingDegreeCriterion ? "(Default)" : " ") << std::endl;
  std::cout<<"Minumum Number of Derived Points (p) : " << derivingMinimumPoint << " " << (!f_derivingMinimumPoint ? "(Default)" : " ") << std::endl;
  std::cout<<"Incident Particle's Velocity (r.SOL) (V) : " << velocity << " " << (!f_velocity ? "(Default)" : " ") << std::endl;
  std::cout<<"Minimum Incident Distance (D) : " << startingMinimumDistance << " " << (!f_startingMinimumDistance ? "(Default)" : " ") << std::endl;

  std::cout<<"Output File's Name (o) : "<< outputFilename << " " << (!f_outputFilename ? "(Default)" : " ") << std::endl;
  std::cout<<"===========================================" << std::endl << std::endl;

  app = new TApplication("app", &argc, argv);

  return Experiment(pi2, pi3, outputFilename, derivingDegreeCriterion, derivingMinimumPoint, velocity, startingMinimumDistance);

}
