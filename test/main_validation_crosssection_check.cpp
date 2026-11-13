////////////////////////////////////////////
//
// RooParticle Test Application
// Physics validation : Rutherford differential cross-section by Monte Carlo
//
// Alpha (4 amu, 2e) at 0.05c on a fixed Au (79e) target. Every event starts at
// the same distance r0 = 1000 fm from the target (so every event has the same
// energy) with velocity (v0, 0, 0) and nominal impact parameter b. Then
//   v_inf = sqrt(v0^2 + 2 kappa/(m r0)),  b_inf = b v0 / v_inf   (exact),
// so a beam uniform over the disc b_inf < b_max (uniform in b^2 -- NOT in b)
// is a beam uniform in b_nom^2.
//
// Each event : RK45 (default tolerance), call step 400 fm/c, stopped by the
// inspectors "VRG" > 0 and "DSTG" > 1000 fm ; theta = angle between the
// asymptote (class asymptote) of the initial and of the final state.
//
// Expected count in a bin [theta1, theta2] (exact, no bin-centre approximation):
//   N (b(theta1)^2 - b(theta2)^2) / b_max^2 ,  b(theta) = a cot(theta/2), a = kappa/(m v_inf^2)
// i.e. the integral of dsigma/dOmega = (a/2)^2 / sin^4(theta/2) over the bin,
// divided by the beam area. Compared by Pearson chi-square and per-bin pulls.
// A second sample drawn uniform in b (the wrong weighting) must be rejected.
//
// No arguments, fixed seed. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"
#include "asymptote.hpp"

#include "iostream"
#include "vector"
#include <cstdio>
#include <ctime>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

static coulombForce gCoulomb;
static Int_t gFail = 0;
static Int_t gCheck = 0;

void Check(const char* what, Double_t val, Double_t tol){
  Bool_t ok = TMath::Finite(val) && val <= tol;
  printf("  %-60s %11.3e < %8.1e  %s\n", what, val, tol, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
void CheckAbove(const char* what, Double_t val, Double_t low){
  Bool_t ok = TMath::Finite(val) && val >= low;
  printf("  %-60s %11.3e > %8.1e  %s\n", what, val, low, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}

static const Double_t kM = 4., kQ = 2., kQT = 79., kV0 = 0.05, kR0 = 1000.;
// RK45 call step : it only caps the adaptive substep, and the inspectors are
// looked at once per call (the overshoot past 1000 fm does not matter, the
// asymptote is taken from wherever the run stops)
static const Double_t kCallDt = 400.;
static Long_t gSub = 0;

// one event : simulated scattering angle (rad) for the nominal impact parameter b
Double_t Simulate(Double_t b, Double_t kap){
  event ev;
  ev.AddForce(&gCoulomb);
  particle* alpha  = ev.AddParticle(new EMparticle(kM, kQ, TVector3(-TMath::Sqrt(kR0*kR0-b*b), b, 0.),
                                                   TVector3(kV0, 0., 0.), false, false));
  particle* target = ev.AddParticle(new EMparticle(197., kQT, TVector3(), TVector3(), true, false));
  ev.SetIntegrator(event::kRK45);
  ev.AddInspector(new inspectorP(alpha, target, "VRG", 0.));
  ev.AddInspector(new inspectorP(alpha, target, "DSTG", kR0));
  asymptote in(alpha, target, kap);
  ev.DeriveInspect(1, kCallDt);
  gSub += ev.GetRK45NSubsteps();
  asymptote out(alpha, target, kap);
  return asymptote::Angle(in.GetVInRel(), out.GetVOutRel());
}

int main(){
  std::clock_t c0 = std::clock();
  const Double_t kap  = gCoulomb.GetK()*kQ*kQT;
  const Double_t vinf = TMath::Sqrt(kV0*kV0 + 2.*kap/(kM*kR0));
  const Double_t a    = kap/(kM*vinf*vinf);
  const Double_t sc   = kV0/vinf;                       // b_inf = sc * b_nom
  const Double_t thMin= 18.*TMath::DegToRad();          // beam edge
  const Double_t bmaxInf = a/TMath::Tan(0.5*thMin);
  const Double_t bmaxNom = bmaxInf/sc;
  const Int_t    N  = 20000;
  const Int_t    NW = 2000;

  printf("validation : Rutherford dsigma/dOmega, alpha 0.05c on fixed Au, Monte Carlo over the beam disc\n");
  printf("  v_inf = %.6f c, a = kappa/(m v_inf^2) = %.4f fm, b_max = %.3f fm (theta_min = 18 deg)\n", vinf, a, bmaxInf);
  printf("  beam area pi b_max^2 = %.2f fm^2 = %.2f b ; N = %d events (uniform in b^2), seed 4357\n",
         TMath::Pi()*bmaxInf*bmaxInf, TMath::Pi()*bmaxInf*bmaxInf/100., N);

  const Double_t edges[] = {20., 25., 30., 40., 50., 60., 75., 90., 110., 135., 180.};
  const Int_t nb = sizeof(edges)/sizeof(edges[0]) - 1;
  std::vector<Double_t> cnt(nb, 0.), cntW(nb, 0.);

  TRandom3 rng(4357);
  Double_t maxEvtErr = 0.;
  Int_t under = 0;
  for(Int_t i=0; i<N; i++){
    Double_t b  = bmaxNom*TMath::Sqrt(rng.Uniform());   // uniform over the disc
    Double_t th = Simulate(b, kap);
    Double_t thR = 2.*TMath::ATan(a/(sc*b));
    maxEvtErr = TMath::Max(maxEvtErr, TMath::Abs(th-thR));
    Double_t d = th*TMath::RadToDeg();
    Int_t k = -1;
    for(Int_t j=0; j<nb; j++) if(d>=edges[j] && d<edges[j+1]) k = j;
    if(k>=0) cnt[k] += 1.; else under++;
  }
  for(Int_t i=0; i<NW; i++){
    Double_t b  = bmaxNom*rng.Uniform();                 // WRONG : uniform in b
    Double_t d  = Simulate(b, kap)*TMath::RadToDeg();
    for(Int_t j=0; j<nb; j++) if(d>=edges[j] && d<edges[j+1]) cntW[j] += 1.;
  }

  printf("\n  %8s %8s %8s %9s %7s | %12s %12s %12s %8s\n", "theta1", "theta2", "N obs", "N exp", "pull",
         "meas mb/sr", "bin-avg mb/sr", "centre mb/sr", "ratio");
  Double_t chi2 = 0., chi2W = 0., maxPull = 0.;
  for(Int_t j=0; j<nb; j++){
    Double_t t1 = edges[j]*TMath::DegToRad(), t2 = edges[j+1]*TMath::DegToRad();
    Double_t b1 = a/TMath::Tan(0.5*t1), b2 = (edges[j+1]>=180.) ? 0. : a/TMath::Tan(0.5*t2);
    Double_t frac = (b1*b1 - b2*b2)/(bmaxInf*bmaxInf);
    Double_t e  = N*frac;
    Double_t pull = (cnt[j]-e)/TMath::Sqrt(e);
    chi2 += pull*pull;
    Double_t eW = NW*frac;
    chi2W += (cntW[j]-eW)*(cntW[j]-eW)/eW;
    maxPull = TMath::Max(maxPull, TMath::Abs(pull));
    Double_t dOmega = 2.*TMath::Pi()*(TMath::Cos(t1)-TMath::Cos(t2));
    Double_t meas = cnt[j]/N*TMath::Pi()*bmaxInf*bmaxInf/dOmega*10.;     // fm^2 -> mb
    Double_t avg  = TMath::Pi()*(b1*b1-b2*b2)/dOmega*10.;
    Double_t sc2  = TMath::Sin(0.25*(t1+t2));
    Double_t cen  = 0.25*a*a/(sc2*sc2*sc2*sc2)*10.;
    printf("  %8.1f %8.1f %8.0f %9.1f %7.2f | %12.4g %12.4g %12.4g %8.4f\n", edges[j], edges[j+1], cnt[j], e, pull,
           meas, avg, cen, meas/avg);
  }
  Double_t prob = TMath::Prob(chi2, nb);
  printf("\n  uniform in b^2 : chi2 = %.2f for %d bins, p = %.3f ; below 20 deg : %d events\n", chi2, nb, prob, under);
  printf("  uniform in b   : chi2 = %.1f for %d bins (wrong weighting, %d events)\n", chi2W, nb, NW);
  printf("  largest per-event |theta_sim - theta_R(b_inf)| = %.2e rad ; %.1f RK45 substeps per event\n",
         maxEvtErr, Double_t(gSub)/(N+NW));

  Check("per-event angle - exact Rutherford (rad)", maxEvtErr, 3e-8);  // measured 4.5e-9 (RK45 default tol)
  CheckAbove("chi-square probability (uniform in b^2)", prob, 1e-3);
  Check("largest |pull| per bin", maxPull, 4.);
  CheckAbove("wrong weighting (uniform in b) is rejected : chi2/ndf", chi2W/nb, 10.);

  printf("\n%d / %d checks passed  (%.1f s)\n", gCheck-gFail, gCheck, Double_t(std::clock()-c0)/CLOCKS_PER_SEC);
  return gFail==0 ? 0 : 1;
}
