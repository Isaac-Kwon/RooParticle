////////////////////////////////////////////
//
// RooParticle Test Application
// Physics validation : two-body scattering on a FREE (recoiling) target
//
// Exact two-body kinematics (non-relativistic). The projectile (alpha) starts at
// (-1000, b, 0) fm with (v0, 0, 0), the target at rest at the origin. Both move.
//
// The textbook lab formulas hold for a target at rest BEFORE the collision, i.e.
// at t -> -inf. A target at rest at a finite distance was moving at -inf (the
// projectile was already pushing it), so the reference lab frame S' here is the
// one in which the target is at rest at -inf. Its velocity v2_in is found by a
// brute-force backward integration (RK45 to 1e9 fm), not from the asymptote class.
// In S' (theory, by hand):
//   relative motion : mu = m1 m2 / M, E = mu v0^2/2 + kappa/r0, v_inf = sqrt(2E/mu),
//                     b_inf = |r0 x v0| / v_inf, theta = 2 atan(kappa/(mu v_inf^2 b_inf))
//   projectile      : tan(theta_1) = sin(theta) / (cos(theta) + m1/m2)
//                     v_1' = v_inf sqrt(m1^2 + m2^2 + 2 m1 m2 cos(theta)) / M
//   target (recoil) : theta_2 = (pi - theta) / 2 ,  v_2' = 2 m1 v_inf sin(theta/2) / M
//   equal masses    : theta_1 + theta_2 = 90 deg
// Simulation : the integrator under test to a stop at 1000 fm, then the class
// asymptote on the final state for the remaining tail, boosted to S'.
// Conservation : total momentum and energy over each run.
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

static const Double_t kM1 = 4., kQ1 = 2.;

event* Build(Double_t m2, Double_t q2, Double_t x0, Double_t b, Double_t v0){
  event* ev = new event();
  ev->AddForce(&gCoulomb);
  ev->AddParticle(new EMparticle(kM1, kQ1, TVector3(x0, b, 0.), TVector3(v0, 0., 0.), false, false));
  ev->AddParticle(new EMparticle(m2,  q2,  TVector3(0., 0., 0.), TVector3(0., 0., 0.), false, false));
  return ev;
}

// fixed dt (or RK45 with that call step) until the pair recedes beyond rstop
void RunToStop(event* ev, Double_t dt, Double_t rstop){
  particle* p1 = ev->getParticle(0);
  particle* p2 = ev->getParticle(1);
  for(Long_t g=0; g<20000000; g++){
    TVector3 dx = p1->GetX()-p2->GetX();
    TVector3 dv = p1->GetV()-p2->GetV();
    if(dx.Dot(dv)>0. && dx.Mag()>=rstop) return;
    ev->DeriveDT(dt);
  }
}

// RK45, tol 1e-13, call step proportional to the distance, to rfar
void FarRun(event* ev, Double_t rfar){
  particle* p1 = ev->getParticle(0);
  particle* p2 = ev->getParticle(1);
  ev->SetIntegrator(event::kRK45);
  ev->SetTolerance(1e-13, 1e-14);
  Double_t vscale = (p1->GetV()-p2->GetV()).Mag();
  for(Long_t g=0; g<10000000; g++){
    TVector3 dx = p1->GetX()-p2->GetX();
    TVector3 dv = p1->GetV()-p2->GetV();
    Double_t r  = dx.Mag();
    if(dx.Dot(dv)>0. && r>=rfar) return;
    ev->DeriveDT(0.02*r/TMath::Max(dv.Mag(), vscale));
  }
}

// Asymptotic lab velocities from a state far away (after FarRun) : the direction
// of the relative velocity there is converged to ~ a b / R^2 (< 1e-15 here), but
// its magnitude still lacks ~ a/R (2e-8 at 1e9 fm). It is set to the exact v_inf
// (energy conservation, by hand) ; the centre-of-mass velocity is exact.
void FarAsymptote(event* ev, Double_t vinf, TVector3 & v1, TVector3 & v2){
  particle* p1 = ev->getParticle(0);
  particle* p2 = ev->getParticle(1);
  Double_t m1 = p1->GetM0(), m2 = p2->GetM0(), M = m1+m2;
  TVector3 vcm  = (p1->GetV()*m1 + p2->GetV()*m2)*(1./M);
  TVector3 vrel = (p1->GetV()-p2->GetV()).Unit()*vinf;
  v1 = vcm + vrel*(m2/M);
  v2 = vcm - vrel*(m1/M);
}

struct integ{ event::integrator id; const char* name; Double_t dt; };
static const integ kIg[] = {
  {event::kVelocityVerlet, "VERLET", 1.},
  {event::kRK4,            "RK4",    2.},
  {event::kRK45,           "RK45",   20.},
  {event::kRK45,           "BRUTE",  0.}};   // RK45 tol 1e-13 to 1e9 fm, no asymptote class
static const Int_t kNI = 4;

int main(){
  std::clock_t c0 = std::clock();
  printf("validation : alpha (m1 = 4, q1 = 2) on a free target, v0 = 0.05c, x0 = -1000 fm\n");
  printf("angles in the frame S' where the target is at rest at t -> -inf (see header)\n");

  struct tcase{ const char* name; Double_t m2, q2; Double_t b[4]; };
  const tcase tc[] = {{"Au (197, 79)", 197., 79., {5., 20., 50., 100.}},
                      {"alpha (4, 2)", 4.,   2.,  {0.5, 1., 2., 5.}}};
  const Double_t x0 = -1000., v0 = 0.05;
  const Double_t d = TMath::RadToDeg();

  Double_t maxA[kNI] = {0}, maxV[kNI] = {0}, maxP[kNI] = {0}, maxE[kNI] = {0}, maxSum[kNI] = {0};
  Double_t maxNaive = 0., maxV2in = 0., maxUin = 0.;

  for(Int_t it=0; it<2; it++){
    const tcase & c = tc[it];
    const Double_t m2 = c.m2, M = kM1+m2, mu = kM1*m2/M, kap = gCoulomb.GetK()*kQ1*c.q2;
    printf("\n== target %s : mu = %.5f amu, kappa = %.6g amu c^2 fm\n", c.name, mu, kap);
    printf("  %6s %8s %9s %9s %9s %9s %10s | %-6s %10s %10s %10s %10s %10s %10s\n",
           "b_nom", "b_inf", "thCM", "th1", "th2", "v1'", "|v2_in|", "integ",
           "dth1 deg", "dth2 deg", "dv1/v1", "dv2/v2", "|dP|/P", "|dE|/E");
    for(Int_t ib=0; ib<4; ib++){
      const Double_t b = c.b[ib];
      // ---- theory by hand
      TVector3 r0(x0, b, 0.), w0(v0, 0., 0.);
      Double_t E    = 0.5*mu*w0.Mag2() + kap/r0.Mag();
      Double_t vinf = TMath::Sqrt(2.*E/mu);
      Double_t binf = r0.Cross(w0).Mag()/vinf;
      Double_t th   = 2.*TMath::ATan(kap/(mu*vinf*vinf*binf));
      Double_t th1  = TMath::ATan2(TMath::Sin(th), TMath::Cos(th) + kM1/m2);
      Double_t th2  = 0.5*(TMath::Pi()-th);
      Double_t v1p  = vinf*TMath::Sqrt(kM1*kM1 + m2*m2 + 2.*kM1*m2*TMath::Cos(th))/M;
      Double_t v2p  = 2.*kM1*vinf*TMath::Sin(0.5*th)/M;

      // ---- incoming asymptotes : backward brute force
      event* evb = Build(m2, c.q2, x0, b, -v0);
      FarRun(evb, 1e9);
      TVector3 v1in, v2in;
      FarAsymptote(evb, vinf, v1in, v2in);
      v1in = -v1in; v2in = -v2in;
      delete evb;
      TVector3 u = v1in - v2in;          // projectile in S' before the collision
      maxV2in = TMath::Max(maxV2in, v2in.Mag()/v0);
      maxUin  = TMath::Max(maxUin, TMath::Abs(u.Mag()-vinf)/vinf);

      // ---- naive : textbook formulas at the nominal v0, b against the raw stop state
      {
        Double_t thN  = 2.*TMath::ATan(kap/(mu*v0*v0*b));
        Double_t th1N = TMath::ATan2(TMath::Sin(thN), TMath::Cos(thN) + kM1/m2);
        event* ev = Build(m2, c.q2, x0, b, v0);
        ev->SetIntegrator(event::kRK4);
        RunToStop(ev, 2., 1000.);
        Double_t th1raw = asymptote::Angle(TVector3(1.,0.,0.), ev->getParticle(0)->GetV());
        maxNaive = TMath::Max(maxNaive, TMath::Abs(th1raw-th1N)*d);
        delete ev;
      }

      for(Int_t k=0; k<kNI; k++){
        event* ev = Build(m2, c.q2, x0, b, v0);
        ev->SetIntegrator(kIg[k].id);
        TVector3 P0 = ev->GetNetMomentum();
        Double_t E0 = ev->GetNetE();
        TVector3 w1, w2;
        if(kIg[k].dt > 0.){
          RunToStop(ev, kIg[k].dt, 1000.);
          asymptote asy(ev->getParticle(0), ev->getParticle(1), kap);
          w1 = asy.GetV1Out() - v2in;
          w2 = asy.GetV2Out() - v2in;
        }else{
          FarRun(ev, 1e9);
          FarAsymptote(ev, vinf, w1, w2);
          w1 -= v2in;
          w2 -= v2in;
        }
        Double_t dP = (ev->GetNetMomentum()-P0).Mag()/P0.Mag();
        Double_t dE = TMath::Abs(ev->GetNetE()-E0)/E0;
        Double_t s1 = asymptote::Angle(u, w1), s2 = asymptote::Angle(u, w2);
        Double_t e1 = (s1-th1)*d, e2 = (s2-th2)*d;
        Double_t ev1 = (w1.Mag()-v1p)/v1p, ev2 = (w2.Mag()-v2p)/v2p;
        if(k==0) printf("  %6.1f %8.4f %9.4f %9.4f %9.4f %9.6f %10.3e |", b, binf, th*d, th1*d, th2*d, v1p, v2in.Mag());
        else     printf("  %6s %8s %9s %9s %9s %9s %10s |", "", "", "", "", "", "", "");
        printf(" %-6s %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n", kIg[k].name, e1, e2, ev1, ev2, dP, dE);
        maxA[k] = TMath::Max(maxA[k], TMath::Max(TMath::Abs(e1), TMath::Abs(e2)));
        maxV[k] = TMath::Max(maxV[k], TMath::Max(TMath::Abs(ev1), TMath::Abs(ev2)));
        maxP[k] = TMath::Max(maxP[k], dP);
        maxE[k] = TMath::Max(maxE[k], dE);
        if(m2==kM1) maxSum[k] = TMath::Max(maxSum[k], TMath::Abs(s1+s2-0.5*TMath::Pi())*d);
        delete ev;
      }
    }
  }

  printf("\n  %-7s %14s %14s %12s %12s %14s\n", "integ", "max|dth| deg", "max|dv|/v", "max|dP|/P", "max|dE|/E", "eq.m. th1+th2-90");
  for(Int_t k=0; k<kNI; k++)
    printf("  %-7s %14.3e %14.3e %12.3e %12.3e %14.3e\n", kIg[k].name, maxA[k], maxV[k], maxP[k], maxE[k], maxSum[k]);
  printf("\n  pitfall : textbook lab angle at nominal (v0, b) vs raw stop direction, worst case %.4f deg\n", maxNaive);
  printf("            target speed at t -> -inf, |v2_in| / v0, largest %.3e\n", maxV2in);

  // tolerances : about 3-7x the largest value measured. Measured (deg / rel. speed / dE/E) :
  // VERLET dt=1 3.1e-4 / 1.2e-5 / 2.5e-13 ; RK4 dt=2 7.9e-8 / 2.9e-9 / 1.4e-10 ;
  // RK45 default tol 1.5e-7 / 3.2e-9 / 1.4e-9 ; BRUTE 2.3e-12 / 4.8e-14 / 2.2e-13.
  // The worst cases are alpha + alpha at b = 0.5-1 fm (closest approach 2.6 fm).
  const Double_t tolA[] = {2e-3, 5e-7, 1e-6, 1e-10};
  const Double_t tolV[] = {5e-5, 2e-8, 2e-8, 1e-12};
  const Double_t tolE[] = {1e-11, 1e-9, 1e-8, 1e-11};
  for(Int_t k=0; k<kNI; k++){
    Check(Form("%-6s lab angles (projectile, recoil) - theory (deg)", kIg[k].name), maxA[k], tolA[k]);
    Check(Form("%-6s lab final speeds - theory (relative)", kIg[k].name), maxV[k], tolV[k]);
    Check(Form("%-6s total momentum drift |dP|/|P|", kIg[k].name), maxP[k], 1e-13);
    Check(Form("%-6s energy drift |dE|/E", kIg[k].name), maxE[k], tolE[k]);
    Check(Form("%-6s equal masses : |theta_1 + theta_2 - 90| (deg)", kIg[k].name), maxSum[k], tolA[k]);
  }
  Check("relative speed before the collision |u| - v_inf (relative)", maxUin, 1e-10);
  CheckAbove("naive textbook-at-nominal error is real (deg)", maxNaive, 0.05);

  printf("\n%d / %d checks passed  (%.1f s)\n", gCheck-gFail, gCheck, Double_t(std::clock()-c0)/CLOCKS_PER_SEC);
  return gFail==0 ? 0 : 1;
}
