////////////////////////////////////////////
//
// RooParticle Test Application
// Physics validation : Coulomb + short-range nuclear force (Woods-Saxon)
//
// Alpha (4 amu, 2e) at 0.12c (E = 26.8 MeV, above the Coulomb barrier) on a
// fixed Au target, starting at r0 = 1000 fm. Forces : coulombForce, and
// coulombForce + woodsSaxonForce (default V0 = 50 MeV, R = 8.886 fm, a = 0.65 fm).
// At 1000 fm the nuclear potential is exp(-1525) ~ 0, so the class asymptote
// (pure 1/r) is exact at the start and at the stop of both runs.
//
// There is no closed form with the nuclear force. Checked :
//   * Coulomb only : exact Rutherford at every b (hand formula at b_inf, v_inf).
//   * Coulomb + WS, large b (Coulomb closest approach far outside R) : the same
//     Rutherford angle, the nuclear force must be invisible.
//   * Coulomb + WS, small b : the angle clearly differs (nuclear attraction).
//   * the b where the deviation sets in (information), against the Coulomb
//     closest approach r_min(b) = a + sqrt(a^2 + b^2).
//   * energy conservation with the combined force (RK45 and VERLET), and the
//     combined-force angle by two different integrators (RK45 vs VERLET, step
//     halving) since there is no formula.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "asymptote.hpp"

#include "iostream"
#include <cstdio>
#include <ctime>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"

static coulombForce    gCoulomb;
static woodsSaxonForce gWS;
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

static const Double_t kM = 4., kQ = 2., kQT = 79., kV0 = 0.12, kR0 = 1000.;

struct runOut{ Double_t theta, edrift, rmin; };

// integ : RK45 (tol 1e-10, call step dt) or a fixed-step integrator with step dt
runOut Run(Double_t b, Bool_t nuclear, event::integrator id, Double_t dt){
  const Double_t kap = gCoulomb.GetK()*kQ*kQT;
  event ev;
  ev.AddForce(&gCoulomb);
  if(nuclear) ev.AddForce(&gWS);
  particle* p1 = ev.AddParticle(new EMparticle(kM, kQ, TVector3(-TMath::Sqrt(kR0*kR0-b*b), b, 0.),
                                               TVector3(kV0, 0., 0.), false, false));
  particle* p2 = ev.AddParticle(new EMparticle(197., kQT, TVector3(), TVector3(), true, false));
  ev.SetIntegrator(id);
  if(id==event::kRK45) ev.SetTolerance(1e-10, 1e-12);
  asymptote in(p1, p2, kap);
  Double_t e0 = ev.GetNetE();
  runOut o; o.edrift = 0.; o.rmin = 1e300;
  for(Long_t g=0; g<50000000; g++){
    TVector3 dx = p1->GetX()-p2->GetX();
    Double_t r = dx.Mag();
    if(r < o.rmin) o.rmin = r;
    if(dx.Dot(p1->GetV())>0. && r>=kR0) break;
    ev.DeriveDT(dt);
    Double_t e = TMath::Abs(ev.GetNetE()-e0)/e0;
    if(e > o.edrift) o.edrift = e;
  }
  asymptote out(p1, p2, kap);
  o.theta = asymptote::Angle(in.GetVInRel(), out.GetVOutRel());
  return o;
}

int main(){
  std::clock_t c0 = std::clock();
  const Double_t kap  = gCoulomb.GetK()*kQ*kQT;
  const Double_t vinf = TMath::Sqrt(kV0*kV0 + 2.*kap/(kM*kR0));
  const Double_t a    = kap/(kM*vinf*vinf);
  const Double_t sc   = kV0/vinf;      // b_inf = sc b
  const Double_t d    = TMath::RadToDeg();
  const Double_t mev  = 931.494;

  printf("validation : alpha 0.12c on fixed Au, Coulomb vs Coulomb + Woods-Saxon\n");
  printf("  E_inf = %.3f MeV, a = %.4f fm, WS : V0 = %.2f MeV, R = %.3f fm, diffuseness %.3f fm\n",
         0.5*kM*vinf*vinf*mev, a, gWS.GetDepth()*mev, gWS.GetRadius(), gWS.GetDiffuseness());
  printf("  Coulomb barrier height at r = R : %.2f MeV (Coulomb + WS : max over r of k q1 q2/r - V0 f(r))\n",
         kap/gWS.GetRadius()*mev);

  const Double_t bs[] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 17., 20., 25., 30., 40., 50., 100.};
  const Int_t nb = sizeof(bs)/sizeof(bs[0]);
  printf("\n  %6s %8s %9s %12s %12s %12s %12s %10s\n", "b_nom", "b_inf", "rmin_C", "theta_R deg",
         "C - R (rad)", "C+WS deg", "C+WS - R deg", "|dE/E|");
  Double_t maxC = 0., maxLarge = 0., minSmall = 1e300, maxE = 0.;
  Double_t onset3 = -1., onset01 = -1., rOn3 = 0., rOn01 = 0.;
  for(Int_t i=0; i<nb; i++){
    Double_t b = bs[i], binf = sc*b;
    Double_t thR  = (b>0.) ? 2.*TMath::ATan(a/binf) : TMath::Pi();
    Double_t rmC  = a + TMath::Sqrt(a*a + binf*binf);
    runOut c  = Run(b, kFALSE, event::kRK45, 100.);
    runOut cw = Run(b, kTRUE,  event::kRK45, 100.);
    Double_t dev = cw.theta - thR;
    printf("  %6.1f %8.4f %9.4f %12.6f %12.2e %12.6f %12.4e %10.2e\n", b, binf, rmC, thR*d, c.theta-thR,
           cw.theta*d, dev*d, cw.edrift);
    maxC = TMath::Max(maxC, TMath::Abs(c.theta-thR));
    maxE = TMath::Max(maxE, TMath::Max(cw.edrift, c.edrift));
    if(b >= 30.) maxLarge = TMath::Max(maxLarge, TMath::Abs(dev));
    // (b = 0 is 180 deg with or without the nuclear force, by symmetry)
    if(b > 0. && b <= 7.) minSmall = TMath::Min(minSmall, TMath::Abs(dev)*d);
    // onset : largest b with a deviation above 1e-3 relative / above 0.1 deg
    if(TMath::Abs(dev) > 1e-3*thR && b > onset3){ onset3 = b; rOn3 = rmC; }
    if(TMath::Abs(dev)*d > 0.1 && b > onset01){ onset01 = b; rOn01 = rmC; }
  }
  printf("\n  deviation sets in : |dtheta| > 1e-3 theta up to b = %.0f fm (Coulomb r_min = %.2f fm = R + %.1f diffuseness)\n",
         onset3, rOn3, (rOn3-gWS.GetRadius())/gWS.GetDiffuseness());
  printf("                      |dtheta| > 0.1 deg   up to b = %.0f fm (Coulomb r_min = %.2f fm = R + %.1f diffuseness)\n",
         onset01, rOn01, (rOn01-gWS.GetRadius())/gWS.GetDiffuseness());

  // no formula with the nuclear force : two integrators, step halving
  printf("\n  Coulomb + WS by VERLET (dt = 0.4, 0.2 fm/c) against RK45 (tol 1e-10) and RK45 tol 1e-12\n");
  printf("  %6s %14s %12s %12s %12s %12s %10s\n", "b_nom", "RK45 deg", "RK45-12 - RK45", "V0.4 - RK45", "V0.2 - RK45",
         "ratio", "V0.2 |dE/E|");
  const Double_t bv[] = {2., 5., 8., 11.};
  Double_t maxV = 0., minRatio = 1e300, maxRatio = 0., maxEV = 0., maxT = 0.;
  for(Int_t i=0; i<4; i++){
    runOut r  = Run(bv[i], kTRUE, event::kRK45, 100.);
    event::integrator id = event::kRK45;
    // RK45 at a tighter tolerance : same function, tolerance through a local event
    runOut rt;
    {
      event ev;
      ev.AddForce(&gCoulomb); ev.AddForce(&gWS);
      particle* p1 = ev.AddParticle(new EMparticle(kM, kQ, TVector3(-TMath::Sqrt(kR0*kR0-bv[i]*bv[i]), bv[i], 0.),
                                                   TVector3(kV0, 0., 0.), false, false));
      particle* p2 = ev.AddParticle(new EMparticle(197., kQT, TVector3(), TVector3(), true, false));
      ev.SetIntegrator(id);
      ev.SetTolerance(1e-12, 1e-14);
      asymptote in(p1, p2, kap);
      for(Long_t g=0; g<50000000; g++){
        TVector3 dx = p1->GetX()-p2->GetX();
        if(dx.Dot(p1->GetV())>0. && dx.Mag()>=kR0) break;
        ev.DeriveDT(100.);
      }
      asymptote out(p1, p2, kap);
      rt.theta = asymptote::Angle(in.GetVInRel(), out.GetVOutRel());
    }
    runOut v4 = Run(bv[i], kTRUE, event::kVelocityVerlet, 0.4);
    runOut v2 = Run(bv[i], kTRUE, event::kVelocityVerlet, 0.2);
    Double_t ratio = (v4.theta-rt.theta)/(v2.theta-rt.theta);
    printf("  %6.1f %14.8f %12.2e %12.2e %12.2e %12.2f %10.2e\n", bv[i], r.theta*d, (rt.theta-r.theta)*d,
           (v4.theta-rt.theta)*d, (v2.theta-rt.theta)*d, ratio, v2.edrift);
    maxT = TMath::Max(maxT, TMath::Abs(rt.theta-r.theta)*d);
    maxV = TMath::Max(maxV, TMath::Abs(v2.theta-rt.theta)*d);
    maxEV = TMath::Max(maxEV, v2.edrift);
    minRatio = TMath::Min(minRatio, ratio);
    maxRatio = TMath::Max(maxRatio, ratio);
  }

  printf("\n");
  Check("Coulomb only : |theta - theta_R| (rad), all b", maxC, 1e-8);
  Check("Coulomb + WS, b >= 30 fm : |theta - theta_R| (rad)", maxLarge, 1e-8);
  CheckAbove("Coulomb + WS, 0 < b <= 7 fm : smallest |theta - theta_R| (deg)", minSmall, 1.);
  Check("energy drift, RK45 (tol 1e-10), both forces", maxE, 1e-8);
  Check("RK45 tol 1e-10 vs 1e-12, Coulomb + WS (deg)", maxT, 1e-5);
  Check("VERLET dt 0.2 vs RK45 tol 1e-12, Coulomb + WS (deg)", maxV, 0.05);
  Check("VERLET dt 0.2 energy drift (max over run), Coulomb + WS", maxEV, 1e-4);  // measured 3e-5
  CheckAbove("VERLET step halving ratio, smallest (2nd order : 4)", minRatio, 3.);
  Check("VERLET step halving ratio, largest (2nd order : 4)", maxRatio, 5.);

  printf("\n%d / %d checks passed  (%.1f s)\n", gCheck-gFail, gCheck, Double_t(std::clock()-c0)/CLOCKS_PER_SEC);
  return gFail==0 ? 0 : 1;
}
