////////////////////////////////////////////
//
// RooParticle tutorial, part 5 : many events
//
// Monte Carlo over the impact parameter for alpha on a fixed gold-like
// target, from a template event and its copies. Writes a TTree and the run
// metadata to Data/<title>_<jobIndex>.root, then reads the file back.
//
// Usage
//   test_tutorial_montecarlo_check                    title "tutorial", job 0
//   test_tutorial_montecarlo_check <title> <jobIndex> (the form macro/runmulti.py uses)
//
// Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"
#include "asymptote.hpp"
#include "runinfo.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TMath.h"
#include "TVector3.h"
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv){
  Int_t nfail = 0;

  //! [args]
  // runmulti.py calls  <program> <title> [extra args] <jobIndex>
  TString title    = (argc > 1) ? argv[1] : "tutorial";
  Int_t   jobIndex = (argc > 2) ? std::atoi(argv[argc - 1]) : 0;
  UInt_t  seed     = runinfo::SeedFromJob(jobIndex);   // 65539 + jobIndex
  //! [args]

  //! [template]
  // Template event: everything the events share. It is never derived itself.
  coulombForce coulomb;                                  // shared by all copies, not owned
  event templ;
  templ.AddForce(&coulomb);
  templ.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  templ.SetIntegrator("VERLET");
  const Double_t dt = 1.;                                // fm/c
  const Double_t kappa = coulomb.GetK() * 2. * 79.;
  //! [template]

  //! [tree]
  gSystem->mkdir("Data", kTRUE);
  TString fname = TString::Format("Data/%s_%d.root", title.Data(), jobIndex);
  TFile file(fname, "RECREATE");
  TTree tree("rutherford", "alpha on fixed Au");

  Double_t imp, sangle, thetaAsy, bEff, vInf, rMin;
  Int_t    nstep, stop;
  tree.Branch("imp",      &imp,      "imp/D");       // nominal impact parameter [fm]
  tree.Branch("SAngle",   &sangle,   "SAngle/D");    // raw angle at the stop [rad]
  tree.Branch("thetaAsy", &thetaAsy, "thetaAsy/D");  // asymptotic angle [rad]
  tree.Branch("bEff",     &bEff,     "bEff/D");      // effective impact parameter [fm]
  tree.Branch("vInf",     &vInf,     "vInf/D");      // speed at infinity [c]
  tree.Branch("rMin",     &rMin,     "rMin/D");      // closest approach [fm]
  tree.Branch("nstep",    &nstep,    "nstep/I");
  tree.Branch("stop",     &stop,     "stop/I");      // event::stopReason
  //! [tree]

  //! [loop]
  const Int_t    nevents = 200;
  const Double_t bmax    = 200.;                     // fm
  TRandom3 rnd(seed);
  Int_t nbad = 0;
  Double_t worstFormula = 0.;
  for(Int_t i = 0; i < nevents; i++){
    // uniform in the area of the beam: b = bmax sqrt(u)
    imp = bmax * TMath::Sqrt(rnd.Uniform());

    event* ev = new event(templ);                    // copies particles, inspectors, settings
    particle* target = ev->getParticle(0);
    particle* alpha  = ev->AddParticle(
        new EMparticle(4., 2., TVector3(-1000., imp, 0.), TVector3(0.05, 0., 0.), false, false));
    ev->AddInspector(new inspectorP(target, alpha, "VRG", 0.));      // receding
    ev->AddInspector(new inspectorP(target, alpha, "DSTG", 1000.));  // back at 1000 fm

    ev->DeriveInspect(1, dt);

    asymptote asy(alpha, target, kappa);
    sangle   = asymptote::Angle(TVector3(1., 0., 0.), alpha->GetV());
    thetaAsy = asy.GetThetaCM();
    bEff     = asy.GetImpactParameter();
    vInf     = asy.GetVinf();
    rMin     = asy.GetClosestApproach();
    nstep    = ev->GetDeriveN();
    stop     = ev->GetStopReason();
    tree.Fill();

    if(stop != event::kStopInspector) nbad++;
    worstFormula = TMath::Max(worstFormula,
        TMath::Abs(thetaAsy - asymptote::RutherfordAngle(bEff, vInf, asy.GetMu(), kappa)));
    delete ev;                                       // deletes its particles and inspectors
  }
  //! [loop]

  //! [meta]
  runinfo::WriteCommon(seed, dt, templ.GetIntegratorName().Data(), nevents, "VRG>0 AND DSTG>1000");
  runinfo::Write("jobIndex", (Long64_t) jobIndex);
  runinfo::Write("projectile", "alpha m=4 q=2 x0=(-1000,b,0) fm v0=(0.05,0,0) c");
  runinfo::Write("target", "fixed m=197 q=79 at origin");
  runinfo::Write("bmax", bmax);
  tree.Write();
  file.Close();
  //! [meta]

  //! [readback]
  TFile in(fname);
  TTree* t = (TTree*) in.Get("rutherford");
  TNamed* s = (TNamed*) in.Get("seed");
  TNamed* g = (TNamed*) in.Get("git_sha");
  std::printf("file %s : %lld events, seed %s, git %s\n", fname.Data(),
              t ? t->GetEntries() : -1LL, s ? s->GetTitle() : "?", g ? g->GetTitle() : "?");
  t->Scan("imp:bEff:SAngle:thetaAsy:rMin:nstep", "", "colsize=11 precision=8", 5);
  //! [readback]

  std::printf("events not ended by the inspectors : %d\n", nbad);
  std::printf("max |thetaAsy - Rutherford(bEff, vInf)| : %.1e rad\n", worstFormula);
  if(nbad != 0) nfail++;
  if(!(worstFormula < 1e-9)) nfail++;
  if(!t || t->GetEntries() != nevents) nfail++;
  if(!s || TString(s->GetTitle()) != TString::Format("%u", seed)) nfail++;
  if(runinfo::SeedFromJob(0) != 65539 || runinfo::SeedFromJob(3) != 65542 || runinfo::SeedFromJob(3, 7) != 7) nfail++;

  std::printf("%s\n", nfail == 0 ? "TUTORIAL 5 OK" : "TUTORIAL 5 FAILED");
  return nfail == 0 ? 0 : 1;
}
