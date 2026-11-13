////////////////////////////////////////////
//
// RooParticle Test Application
// Physics validation : two-body scattering on a fixed Coulomb centre
//
// Exact two-body theory (non-relativistic unless stated), every theory value
// computed BY HAND in this file from the initial energy and angular momentum
// (never from the nominal b, v0 of the set-up) :
//
//   E = m v0^2/2 + kappa/r0 ,  v_inf = sqrt(2E/m) ,  b_inf = |r0 x v0| / v_inf
//   a = kappa / (m v_inf^2) ,  theta_R = 2 atan(a / b_inf)
//   r_min = a + sqrt(a^2 + b_inf^2)   (head-on : r_min = 2a = kappa/E)
//
// 0. pitfall      : nominal vs asymptotic b, v and Rutherford angle at a
//                   finite start / stop distance (information + checks).
// 1. angle vs b   : every integrator, b = 0 ... 500 fm, two energies.
//                   Simulated angle = between the asymptote of the initial
//                   state and that of the final state (class asymptote).
// 1b. independent : brute-force RK45 out to 1e9 fm in both time directions
//                   (no asymptote class) against the hand formula.
// 2. closest      : head-on and off-axis distance of closest approach from
//    approach       the trajectory minimum and from path::GetDCA (discrete
//                   and interpolated, fine and coarse record), and the
//                   path-to-path GetDCA of a free pair.
// 3. relativistic : SetRelativistic(kTRUE) against the exact relativistic
//                   fixed-centre deflection (Darwin) and, as information,
//                   against the non-relativistic Rutherford formula.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "path.hpp"
#include "asymptote.hpp"

#include "iostream"
#include "vector"
#include <cstdio>
#include <ctime>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"

static coulombForce gCoulomb;
static Int_t gFail = 0;
static Int_t gCheck = 0;

static const Double_t kMAlpha = 4.;      // amu
static const Double_t kQAlpha = 2.;      // e
static const Double_t kMAu    = 197.;
static const Double_t kQAu    = 79.;
static const Double_t kMPositron = 5.48580e-4;

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

Double_t Kappa(Double_t q1, Double_t q2){ return gCoulomb.GetK()*q1*q2; }

// ---- exact two-body theory, by hand from one state (r, v) of the relative motion
//! [hand_theory]
struct theory{
  Double_t E, vinf, L, binf, a, theta, rmin;
};
theory Hand(TVector3 r, TVector3 v, Double_t mu, Double_t kappa){
  theory t;
  t.E    = 0.5*mu*v.Mag2() + kappa/r.Mag();
  t.vinf = TMath::Sqrt(2.*t.E/mu);
  t.L    = mu*r.Cross(v).Mag();
  t.binf = t.L/(mu*t.vinf);
  t.a    = kappa/(mu*t.vinf*t.vinf);
  t.theta= (t.binf>0.) ? 2.*TMath::ATan(t.a/t.binf) : TMath::Pi();
  t.rmin = t.a + TMath::Sqrt(t.a*t.a + t.binf*t.binf);
  return t;
}
//! [hand_theory]

// Exact relativistic deflection by a fixed centre, H = sqrt(p^2+m^2) + kappa/r
// (c = 1). Orbit u'' + G^2 u = -kappa W / L^2, G^2 = 1 - kappa^2/L^2 :
//   theta = pi - (2/G) atan(G p_inf L / (kappa W)),  p_inf = sqrt(W^2 - m^2).
//! [darwin]
Double_t DarwinAngle(Double_t W, Double_t m, Double_t L, Double_t kappa){
  Double_t pinf = TMath::Sqrt(W*W - m*m);
  Double_t G    = TMath::Sqrt(1. - kappa*kappa/(L*L));
  return TMath::Pi() - (2./G)*TMath::ATan(G*pinf*L/(kappa*W));
}
//! [darwin]

// projectile = particle 0, fixed target = particle 1
event* Build(Double_t m1, Double_t q1, Double_t x0, Double_t b, Double_t v0,
             Bool_t record=kFALSE, Bool_t fixed2=kTRUE, Double_t m2=kMAu, Double_t q2=kQAu){
  event* ev = new event();
  ev->AddForce(&gCoulomb);
  ev->AddParticle(new EMparticle(m1, q1, TVector3(x0, b, 0.), TVector3(v0, 0., 0.), false, record));
  ev->AddParticle(new EMparticle(m2, q2, TVector3(0., 0., 0.), TVector3(0., 0., 0.), fixed2, record));
  return ev;
}

// Fixed dt until the pair recedes and is farther than rstop.
// rmin  : parabola through the three samples around the smallest distance
// rdisc : smallest sampled distance
struct runInfo{ Double_t rmin = 1e300, rdisc = 1e300; Long_t n = 0; Double_t edrift = 0.; };
runInfo RunToStop(event* ev, Double_t dt, Double_t rstop){
  runInfo ri;
  particle* p1 = ev->getParticle(0);
  particle* p2 = ev->getParticle(1);
  Double_t e0 = ev->GetNetE();
  Double_t rr[3] = {0., 0., 0.};
  for(Long_t g=0; g<20000000; g++){
    TVector3 dx = p1->GetX()-p2->GetX();
    TVector3 dv = p1->GetV()-p2->GetV();
    Double_t r  = dx.Mag();
    rr[0] = rr[1]; rr[1] = rr[2]; rr[2] = r;
    if(r < ri.rdisc) ri.rdisc = r;
    if(g>=2 && rr[1]<=rr[0] && rr[1]<=rr[2]){
      Double_t den = rr[0]-2.*rr[1]+rr[2];
      Double_t rv  = (den>0.) ? rr[1]-(rr[2]-rr[0])*(rr[2]-rr[0])/(8.*den) : rr[1];
      if(rv < ri.rmin) ri.rmin = rv;
    }
    if(dx.Dot(dv)>0. && r>=rstop) break;
    ev->DeriveDT(dt);
    ri.n++;
  }
  if(ri.rdisc < ri.rmin) ri.rmin = ri.rdisc;
  ri.edrift = TMath::Abs(ev->GetNetE()-e0)/TMath::Abs(e0);
  return ri;
}

// RK45 out to distance rfar (call step proportional to the distance).
// Integrates the relativistic dynamics if the event is set so.
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

// brute-force angle between the incoming and outgoing asymptotes of the
// projectile, from the state (x0, b, v0) : backward run (reversed velocity)
// and forward run, both to rfar. No asymptote class involved.
Double_t BruteAngle(Double_t m1, Double_t q1, Double_t x0, Double_t b, Double_t v0,
                    Double_t rfar, Bool_t rel, Double_t* speedOut=nullptr, Double_t* edrift=nullptr,
                    Double_t* rOut=nullptr){
  event* evb = Build(m1, q1, x0, b, -v0);
  event* evf = Build(m1, q1, x0, b,  v0);
  evb->SetRelativistic(rel);
  evf->SetRelativistic(rel);
  Double_t e0 = rel ? evf->GetNetERel() : evf->GetNetE();
  FarRun(evb, rfar);
  FarRun(evf, rfar);
  TVector3 vin  = -evb->getParticle(0)->GetV();
  TVector3 vout =  evf->getParticle(0)->GetV();
  if(speedOut) *speedOut = vout.Mag();
  if(rOut) *rOut = evf->getParticle(0)->GetX().Mag();
  if(edrift){
    Double_t e1 = rel ? evf->GetNetERel() : evf->GetNetE();
    *edrift = TMath::Abs(e1-e0)/TMath::Abs(e0);
  }
  delete evb; delete evf;
  return asymptote::Angle(vin, vout);
}

// ======================================================================
void Pitfall(){
  printf("\n== 0. finite-distance pitfall : alpha 0.05c on fixed Au, nominal b = 50 fm, x0 = -1000 fm\n");
  const Double_t kap = Kappa(kQAlpha, kQAu);
  const Double_t x0 = -1000., b = 50., v0 = 0.05;
  theory th = Hand(TVector3(x0, b, 0.), TVector3(v0, 0., 0.), kMAlpha, kap);
  Double_t thNom = 2.*TMath::ATan(kap/(kMAlpha*v0*v0*b));

  event* ev = Build(kMAlpha, kQAlpha, x0, b, v0);
  ev->SetIntegrator(event::kRK4);
  asymptote asyIn(ev->getParticle(0), ev->getParticle(1), kap);
  RunToStop(ev, 1., 1000.);
  asymptote asyOut(ev->getParticle(0), ev->getParticle(1), kap);
  Double_t thRaw = asymptote::Angle(TVector3(v0, 0., 0.), ev->getParticle(0)->GetV());
  Double_t thSim = asymptote::Angle(asyIn.GetVInRel(), asyOut.GetVOutRel());
  delete ev;

  const Double_t d = TMath::RadToDeg();
  printf("  %-46s %14s %14s\n", "", "nominal", "asymptotic");
  printf("  %-46s %14.6f %14.6f\n", "impact parameter b (fm)", b, th.binf);
  printf("  %-46s %14.6f %14.6f\n", "speed v (c)", v0, th.vinf);
  printf("  %-46s %14.6f %14.6f\n", "Rutherford angle at these b, v (deg)", thNom*d, th.theta*d);
  printf("  %-46s %14.6f\n", "raw angle v(start) -> v(stop at 1000 fm) (deg)", thRaw*d);
  printf("  %-46s %14.6f\n", "asymptote(start) -> asymptote(stop) (deg)", thSim*d);
  CheckAbove("nominal formula is off by (deg)", TMath::Abs(thNom-th.theta)*d, 0.5);
  // (smaller than the nominal error : the tilt of the start and the residual
  //  deflection after the stop partly cancel)
  CheckAbove("raw stop direction is off by (deg)", TMath::Abs(thRaw-th.theta)*d, 0.03);
  Check("asymptotic sim - asymptotic theory (deg)", TMath::Abs(thSim-th.theta)*d, 1e-5);
}

// ======================================================================
void AngleVsB(){
  printf("\n== 1. Rutherford angle vs impact parameter, alpha on fixed Au, x0 = -1000 fm, stop at 1000 fm\n");
  printf("   theory : hand formula at (b_inf, v_inf) ; sim : asymptote(initial) -> asymptote(final)\n");
  const Double_t kap = Kappa(kQAlpha, kQAu);
  const Double_t bs[]  = {0., 1., 5., 10., 20., 50., 100., 200., 500.};
  const Int_t    nb    = sizeof(bs)/sizeof(bs[0]);
  const Double_t vs[]  = {0.05, 0.1};
  struct integ{ event::integrator id; const char* name; Double_t dt; Double_t tolDeg; Double_t tolE; };
  // tolerances : about 3-5x the largest error measured over this table
  // (measured : TAYLOR 0.147 deg, 9.0e-3 ; SEULER 4.4e-4, 4.7e-6 ; VERLET 4.1e-4, 9.8e-12 ;
  //  RK4 4.3e-9, 4.0e-12 ; RK45 4.7e-8, 1.0e-9 -- the worst case is always v0 = 0.1c, b ~ 5 fm)
  const integ ig[] = {
    {event::kTaylor,          "TAYLOR", 2.,  0.5,   3e-2},
    {event::kSymplecticEuler, "SEULER", 2.,  2e-3,  2e-5},
    {event::kVelocityVerlet,  "VERLET", 2.,  2e-3,  5e-11},
    {event::kRK4,             "RK4",    2.,  2e-8,  2e-11},
    {event::kRK45,            "RK45",   20., 2e-7,  5e-9}};
  const Int_t ni = sizeof(ig)/sizeof(ig[0]);

  printf("  %5s %7s %8s %8s %12s |", "v0", "b_nom", "b_inf", "v_inf", "theta_R");
  for(Int_t k=0; k<ni; k++) printf(" %10s", ig[k].name);
  printf("   (sim - theory, deg)\n");

  Double_t maxErr[8] = {0}, maxE[8] = {0};
  for(Int_t iv=0; iv<2; iv++){
    for(Int_t ib=0; ib<nb; ib++){
      const Double_t x0 = -1000., b = bs[ib], v0 = vs[iv];
      theory th = Hand(TVector3(x0, b, 0.), TVector3(v0, 0., 0.), kMAlpha, kap);
      printf("  %5.2f %7.1f %8.3f %8.5f %12.6f |", v0, b, th.binf, th.vinf, th.theta*TMath::RadToDeg());
      for(Int_t k=0; k<ni; k++){
        event* ev = Build(kMAlpha, kQAlpha, x0, b, v0);
        ev->SetIntegrator(ig[k].id);
        asymptote asyIn(ev->getParticle(0), ev->getParticle(1), kap);
        runInfo ri = RunToStop(ev, ig[k].dt, 1000.);
        asymptote asyOut(ev->getParticle(0), ev->getParticle(1), kap);
        Double_t th_s = asymptote::Angle(asyIn.GetVInRel(), asyOut.GetVOutRel());
        Double_t err  = (th_s - th.theta)*TMath::RadToDeg();
        printf(" %10.2e", err);
        if(TMath::Abs(err) > maxErr[k] || !TMath::Finite(err)) maxErr[k] = TMath::Finite(err) ? TMath::Abs(err) : 1e300;
        if(ri.edrift > maxE[k]) maxE[k] = ri.edrift;
        delete ev;
      }
      printf("\n");
    }
  }
  printf("\n  %-8s %6s %16s %16s\n", "integ", "dt", "max|dtheta| deg", "max |dE/E|");
  for(Int_t k=0; k<ni; k++) printf("  %-8s %6.1f %16.3e %16.3e\n", ig[k].name, ig[k].dt, maxErr[k], maxE[k]);
  for(Int_t k=0; k<ni; k++){
    Check(Form("%-6s max |theta_sim - theta_R| (deg), 18 cases", ig[k].name), maxErr[k], ig[k].tolDeg);
    Check(Form("%-6s max energy drift |dE/E|", ig[k].name), maxE[k], ig[k].tolE);
  }
}

// ======================================================================
void Independent(){
  printf("\n== 1b. independent of the asymptote class : RK45 (tol 1e-13) to 1e9 fm, backward and forward\n");
  const Double_t kap = Kappa(kQAlpha, kQAu);
  const Double_t bs[] = {5., 50., 500.};
  printf("  %6s %14s %14s %11s %11s %11s %11s\n", "b_nom", "theta_R deg", "brute deg", "d(rad)",
         "(v-v_inf)/v", "(v-v(R))/v", "asy-hand");
  Double_t maxD = 0., maxV = 0., maxA = 0.;
  for(Int_t ib=0; ib<3; ib++){
    const Double_t x0 = -1000., b = bs[ib], v0 = 0.05;
    theory th = Hand(TVector3(x0, b, 0.), TVector3(v0, 0., 0.), kMAlpha, kap);
    Double_t vout, rend;
    Double_t thB = BruteAngle(kMAlpha, kQAlpha, x0, b, v0, 1e9, kFALSE, &vout, nullptr, &rend);
    // the asymptote class from the initial state, against the hand formula
    asymptote asy(TVector3(x0, b, 0.), TVector3(v0, 0., 0.), kMAlpha, kap);
    Double_t dA = TMath::Abs(asy.GetThetaCM()-th.theta);
    Double_t dB = TMath::Abs(thB-th.theta);
    // at the finite end distance R the speed is still v(R) = sqrt(2(E - kappa/R)/m),
    // which differs from v_inf by ~ a/R (2.3e-8 at R = 1e9 fm)
    Double_t vR = TMath::Sqrt(2.*(th.E - kap/rend)/kMAlpha);
    Double_t dVinf = (vout-th.vinf)/th.vinf;
    Double_t dV = TMath::Abs(vout-vR)/vR;
    printf("  %6.1f %14.9f %14.9f %11.2e %11.2e %11.2e %11.2e\n", b, th.theta*TMath::RadToDeg(), thB*TMath::RadToDeg(),
           dB, dVinf, dV, dA);
    maxD = TMath::Max(maxD, dB); maxV = TMath::Max(maxV, dV); maxA = TMath::Max(maxA, dA);
  }
  // residual deflection beyond 1e9 fm is ~ a b / (2 R^2) < 1e-13 rad ; RK45 error ~ 1e-12
  Check("brute-force angle - hand Rutherford (rad)", maxD, 1e-10);
  Check("brute-force speed at R - hand v(R) (relative)", maxV, 1e-10);
  Check("asymptote::GetThetaCM(initial state) - hand (rad)", maxA, 1e-12);
}

// ======================================================================
void ClosestApproach(){
  printf("\n== 2. distance of closest approach, alpha on fixed Au (VERLET and RK4, dt = 1 fm/c)\n");
  printf("   d = a + sqrt(a^2 + b_inf^2) (head-on : kappa/E) ; 'traj' = parabola-refined minimum of r(t)\n");
  printf("   of the run ; GetDCA on the recorded path (every step) and on every 25th sample of it\n");
  const Double_t kap = Kappa(kQAlpha, kQAu);
  printf("  %-6s %5s %5s %10s %10s %10s | %10s %10s %10s %10s\n", "integ", "v0", "b", "d (fm)", "d_nominal",
         "traj-d", "intp1-traj", "disc25-trj", "bound25", "intp25-trj");
  const Double_t vs[] = {0.05, 0.1};
  const Double_t bs[] = {0., 20., 50.};
  const event::integrator ids[] = {event::kVelocityVerlet, event::kRK4};
  const char* names[] = {"VERLET", "RK4"};
  const Double_t tolTraj[] = {3e-5, 1e-10};
  Double_t maxTraj[2] = {0., 0.};
  Double_t maxFine = 0., maxCoarseRel = 0., minDCA = 1e300, maxBoundViol = -1e300;
  for(Int_t k=0; k<2; k++) for(Int_t iv=0; iv<2; iv++) for(Int_t ib=0; ib<3; ib++){
    const Double_t x0 = -1000., b = bs[ib], v0 = vs[iv];
    theory th = Hand(TVector3(x0, b, 0.), TVector3(v0, 0., 0.), kMAlpha, kap);
    Double_t dNom = (b==0.) ? kap/(0.5*kMAlpha*v0*v0) : 0.;
    event* ev = Build(kMAlpha, kQAlpha, x0, b, v0, kTRUE);
    ev->SetIntegrator(ids[k]);
    runInfo ri = RunToStop(ev, 1., 1000.);
    path* pa = ev->getParticle(0)->GetPath();
    // coarse record : every 25th sample of the fine one (25 fm/c apart)
    path coarse;
    for(Int_t i=0; i<pa->GetRecordSize(); i+=25) coarse.AddRecord(pa->GetRecordX(i), pa->GetRecordV(i), pa->GetRecordT(i));
    Double_t dFi, dCd, dCi, t;
    pa->GetDCA(TVector3(0.,0.,0.), dFi, t, false, true);
    coarse.GetDCA(TVector3(0.,0.,0.), dCd, t, false, false);
    coarse.GetDCA(TVector3(0.,0.,0.), dCi, t, false, true);
    // discrete sampling error bound : r(t) ~ rmin + (rdd/2) (t-t0)^2 and a sample
    // lies within T/2 = 12.5 fm/c of t0 ; rdd = v_t^2/rmin + kappa/(m rmin^2), v_t = L/(m rmin)
    Double_t vt    = th.L/(kMAlpha*th.rmin);
    Double_t rdd   = vt*vt/th.rmin + kap/(kMAlpha*th.rmin*th.rmin);
    Double_t bound = 0.5*rdd*12.5*12.5;
    printf("  %-6s %5.2f %5.1f %10.5f %10s %10.2e | %10.2e %10.2e %10.2e %10.2e\n", names[k], v0, b, th.rmin,
           (b==0.) ? Form("%10.5f", dNom) : "-", ri.rmin-th.rmin, dFi-ri.rmin, dCd-ri.rmin, bound, dCi-ri.rmin);
    maxTraj[k]   = TMath::Max(maxTraj[k], TMath::Abs(ri.rmin-th.rmin));
    maxFine      = TMath::Max(maxFine, TMath::Abs(dFi-ri.rmin));
    maxCoarseRel = TMath::Max(maxCoarseRel, TMath::Abs(dCi-ri.rmin)/bound);
    minDCA       = TMath::Min(minDCA, TMath::Min(dFi, dCi));
    // discrete DCA must be above the true minimum and within the bound
    maxBoundViol = TMath::Max(maxBoundViol, TMath::Max(ri.rmin-dCd-1e-9, dCd-ri.rmin-bound));
    delete ev;
  }
  for(Int_t k=0; k<2; k++) Check(Form("%-6s trajectory minimum - d (fm)", names[k]), maxTraj[k], tolTraj[k]);
  Check("GetDCA interpolated, every step - trajectory minimum (fm)", maxFine, 1e-8);
  Check("GetDCA discrete, every 25th : outside [traj, traj+bound] (fm)", TMath::Max(maxBoundViol, 0.), 0.);
  Check("GetDCA interpolated, every 25th : |DCA-traj| / bound", maxCoarseRel, 1e-3);
  CheckAbove("GetDCA is not zero (fixed bug) : smallest DCA (fm)", minDCA, 10.);

  // free pairs : path-to-path DCA, both particles recorded on the same time grid
  printf("\n  free pairs (both bodies move), RK4 dt = 1, path-to-path GetDCA ; d = kappa/E_rel (head-on)\n");
  printf("  %-12s %5s %12s %12s %12s %12s %12s\n", "pair", "v0", "d (fm)", "d_nominal", "traj-d", "intp-traj", "disc-traj");
  struct fcase{ const char* name; Double_t m2, q2, v0; };
  const fcase fc[] = {{"alpha+Au", kMAu, kQAu, 0.05}, {"alpha+Au", kMAu, kQAu, 0.1}, {"alpha+alpha", kMAlpha, kQAlpha, 0.02}};
  Double_t maxPPt = 0., maxPPi = 0.;
  for(Int_t i=0; i<3; i++){
    const Double_t x0 = -1000., v0 = fc[i].v0;
    const Double_t kap2 = Kappa(kQAlpha, fc[i].q2), mu = kMAlpha*fc[i].m2/(kMAlpha+fc[i].m2);
    theory th = Hand(TVector3(x0, 0., 0.), TVector3(v0, 0., 0.), mu, kap2);
    event* ev = Build(kMAlpha, kQAlpha, x0, 0., v0, kTRUE, kFALSE, fc[i].m2, fc[i].q2);
    ev->SetIntegrator(event::kRK4);
    runInfo ri = RunToStop(ev, 1., 1000.);
    Double_t dd, di, t;
    ev->getParticle(0)->GetPath()->GetDCA(ev->getParticle(1)->GetPath(), dd, t, false, false);
    ev->getParticle(0)->GetPath()->GetDCA(ev->getParticle(1)->GetPath(), di, t, false, true);
    printf("  %-12s %5.2f %12.6f %12.6f %12.2e %12.2e %12.2e\n", fc[i].name, v0, th.rmin, kap2/(0.5*mu*v0*v0),
           ri.rmin-th.rmin, di-ri.rmin, dd-ri.rmin);
    maxPPt = TMath::Max(maxPPt, TMath::Abs(ri.rmin-th.rmin));
    maxPPi = TMath::Max(maxPPi, TMath::Abs(di-ri.rmin));
    delete ev;
  }
  Check("free pair : trajectory minimum - kappa/E_rel (fm)", maxPPt, 1e-8);
  Check("free pair : path-to-path GetDCA interpolated - traj. min (fm)", maxPPi, 1e-8);
}

// ======================================================================
void Relativistic(){
  printf("\n== 3. relativistic mode (RK45, tol 1e-13, brute force to 1e9 x scale) on a fixed Au centre\n");
  printf("   exact : Darwin  theta = pi - (2/G) atan(G p L/(kappa W)), G = sqrt(1-kappa^2/L^2)\n");
  printf("   NR    : 2 atan(kappa/(m0 v^2 b)) ;  pv : 2 atan(kappa/(p v b))  (v, p, b at infinity)\n");
  struct rcase{ const char* name; Double_t m, q, v0, thetaPick; };
  const rcase cs[] = {
    {"alpha",    kMAlpha,    kQAlpha, 0.01, 60.},
    {"alpha",    kMAlpha,    kQAlpha, 0.05, 60.},
    {"alpha",    kMAlpha,    kQAlpha, 0.05, 150.},
    {"alpha",    kMAlpha,    kQAlpha, 0.30, 60.},
    {"positron", kMPositron, 1.,      0.50, 60.},
    {"positron", kMPositron, 1.,      0.90, 20.},
    {"positron", kMPositron, 1.,      0.90, 60.},
    {"positron", kMPositron, 1.,      0.99, 60.}};
  const Int_t nc = sizeof(cs)/sizeof(cs[0]);
  printf("  %-8s %5s %9s %6s %11s %11s %10s %10s %10s %10s\n", "proj", "v0", "b (fm)", "k/L",
         "Darwin deg", "sim deg", "sim-D rad", "NR-sim deg", "pv-sim deg", "|dE/E|");
  Double_t maxD = 0., maxE = 0., maxLowRel = 0.;
  const Double_t d = TMath::RadToDeg();
  for(Int_t i=0; i<nc; i++){
    const rcase & c = cs[i];
    const Double_t kap = Kappa(c.q, kQAu);
    Double_t g0 = 1./TMath::Sqrt(1.-c.v0*c.v0);
    // b chosen so that the "pv" formula at the nominal state gives about thetaPick
    Double_t b  = kap/(g0*c.m*c.v0*c.v0*TMath::Tan(0.5*c.thetaPick*TMath::DegToRad()));
    Double_t x0 = -20.*TMath::Max(b, kap/(g0*c.m*c.v0*c.v0));
    // conserved quantities, by hand
    TVector3 r0(x0, b, 0.), p0(g0*c.m*c.v0, 0., 0.);
    Double_t W   = g0*c.m + kap/r0.Mag();
    Double_t L   = r0.Cross(p0).Mag();
    Double_t pin = TMath::Sqrt(W*W - c.m*c.m);
    Double_t vin = pin/W;
    Double_t bin = L/pin;
    Double_t thD = DarwinAngle(W, c.m, L, kap);
    Double_t thNR = 2.*TMath::ATan(kap/(c.m*vin*vin*bin));
    Double_t thPV = 2.*TMath::ATan(kap/(pin*vin*bin));
    Double_t ed;
    Double_t thS = BruteAngle(c.m, c.q, x0, b, c.v0, 1e9*TMath::Max(1., b/1000.), kTRUE, nullptr, &ed);
    printf("  %-8s %5.2f %9.3f %6.3f %11.6f %11.6f %10.2e %10.2e %10.2e %10.2e\n", c.name, c.v0, b, kap/L,
           thD*d, thS*d, thS-thD, (thNR-thS)*d, (thPV-thS)*d, ed);
    maxD = TMath::Max(maxD, TMath::Abs(thS-thD));
    maxE = TMath::Max(maxE, ed);
    if(c.v0 <= 0.05){
      // low speed : the relativistic correction is O(v^2)
      maxLowRel = TMath::Max(maxLowRel, TMath::Abs(thNR-thS)/thS/(vin*vin));
    }
  }
  Check("relativistic sim - exact Darwin angle (rad), all speeds", maxD, 1e-7);
  Check("relativistic energy drift (GetNetERel)", maxE, 1e-9);
  Check("v <= 0.05 : |theta_NR - theta_sim| / (theta v^2)", maxLowRel, 1.);
}

int main(){
  std::clock_t c0 = std::clock();
  printf("validation : two-body scattering on a fixed Coulomb centre (alpha on Au unless stated)\n");
  printf("kappa = k q1 q2, k = %.8g amu c^2 fm / e^2\n", gCoulomb.GetK());
  Pitfall();
  AngleVsB();
  Independent();
  ClosestApproach();
  Relativistic();
  printf("\n%d / %d checks passed  (%.1f s)\n", gCheck-gFail, gCheck, Double_t(std::clock()-c0)/CLOCKS_PER_SEC);
  return gFail==0 ? 0 : 1;
}
