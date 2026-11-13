////////////////////////////////////////////
//
// RooParticle Test Application
// Physics validation : N-body (three or more bodies), Coulomb only
//
// There is no closed form for N >= 3. The evidence here is numerical :
//
// A. conservation and convergence, three configurations
//      C1  3 free bodies (alpha, alpha, 12C)
//      C2  fixed Au target + 2 alphas that also repel each other
//      C3  4 free bodies (p, alpha, 12C, 16O)
//    For every integrator and step : largest energy error over the run,
//    final momentum and angular momentum drift, final position error against
//    a reference (RK4, dt = 0.25 fm/c, cross-checked with RK45 tol 1e-13),
//    and the order of convergence from step halving. (C2 : momentum is not
//    conserved, the fixed target takes it ; angular momentum about the target is.)
//    RK45 is run with a large call step (200 fm/c) so that its tolerance, not
//    the call step, sets the substep.
// B. decoupling : alpha on fixed Au (b = 50 fm, 0.05c) plus a proton at rest
//    at distance D. The pair must reproduce the exact two-body Rutherford angle
//    within the first-order estimate of the perturbation by the (nearly uniform)
//    field F = k q_a q_p / D^2 of the proton over the run time T :
//      est = F/(m v0) [ T + |dtheta/db| Int |x(t)| dt ]
//    (direct impulse, plus the torque about the target that changes L and so b).
//    The deviation must stay below est, be of its size, and scale as 1/D^2.
// C. asymptotic states : RK45 (tol 1e-13) out to 1e12 fm. Outgoing directions
//    at stopping distance R (smallest pair distance) against R = 1e12 fm ;
//    convergence rate k in  delta(R) ~ C / R^k (fit 1e4 ... 1e8 fm).
//      two bodies, fixed centre (or relative direction) : k = 2
//      three or more bodies (also with a fixed target),
//      and two free bodies seen in the lab               : k = 1
//    Speeds converge as 1/R in every case (the potential energy left ~ 1/R).
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
#include "vector"
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
  printf("  %-62s %11.3e < %8.1e  %s\n", what, val, tol, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
void CheckRange(const char* what, Double_t val, Double_t lo, Double_t hi){
  Bool_t ok = TMath::Finite(val) && val >= lo && val <= hi;
  printf("  %-62s %11.4f in [%5.2f, %5.2f]  %s\n", what, val, lo, hi, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}

struct body{ Double_t m, q; TVector3 x, v; Bool_t fixed; };
struct config{ const char* name; std::vector<body> b; Double_t T; Bool_t pCons; };

event* Make(const config & c){
  event* ev = new event();
  ev->AddForce(&gCoulomb);
  for(size_t i=0; i<c.b.size(); i++)
    ev->AddParticle(new EMparticle(c.b[i].m, c.b[i].q, c.b[i].x, c.b[i].v, c.b[i].fixed, false));
  return ev;
}

std::vector<config> Configs(){
  std::vector<config> cs;
  //! [configs]
  config c1; c1.name = "C1 3 free"; c1.T = 20000.; c1.pCons = kTRUE;
  c1.b.push_back({4.,  2., TVector3(-400.,    0.,  0.), TVector3( 0.04, 0.,   0.), kFALSE});
  c1.b.push_back({4.,  2., TVector3( 400.,   20.,  0.), TVector3(-0.04, 0.,   0.), kFALSE});
  c1.b.push_back({12., 6., TVector3(   0., -300., 30.), TVector3( 0.,   0.03, 0.), kFALSE});
  cs.push_back(c1);
  config c2; c2.name = "C2 Au+2a"; c2.T = 40000.; c2.pCons = kFALSE;
  c2.b.push_back({197., 79., TVector3(    0.,     0.,  0.), TVector3(0.,   0.,   0.), kTRUE});
  c2.b.push_back({4.,    2., TVector3(-1000.,    30.,  0.), TVector3(0.05, 0.,   0.), kFALSE});
  c2.b.push_back({4.,    2., TVector3(   20., -1000., 10.), TVector3(0.,   0.05, 0.), kFALSE});
  cs.push_back(c2);
  config c3; c3.name = "C3 4 free"; c3.T = 40000.; c3.pCons = kTRUE;
  c3.b.push_back({1.00728, 1., TVector3(-600.,   10.,   0.), TVector3( 0.03,   0.,    0.), kFALSE});
  c3.b.push_back({4.,      2., TVector3( 500.,  -15.,   5.), TVector3(-0.025,  0.,    0.), kFALSE});
  c3.b.push_back({12.,     6., TVector3(  15., -450., -10.), TVector3( 0.,     0.022, 0.), kFALSE});
  c3.b.push_back({16.,     8., TVector3( -10.,  400.,  20.), TVector3( 0.,    -0.02,  0.), kFALSE});
  cs.push_back(c3);
  //! [configs]
  return cs;
}

// scale of momentum and angular momentum (about the origin) for relative drifts
void Scales(event* ev, Double_t & ps, Double_t & ls){
  ps = 0.; ls = 0.;
  for(Int_t i=0; i<ev->getNParticle(); i++){
    particle* p = ev->getParticle(i);
    if(p->IsInvincible()) continue;
    TVector3 mom = p->GetV()*p->GetM0();
    ps += mom.Mag();
    ls += p->GetX().Cross(mom).Mag();
  }
}

struct result{
  Double_t emax = 0., dp = 0., dl = 0., dx = 0., rmin = 1e300;
  Long_t   nsub = 0;
  std::vector<TVector3> x;
};

// Fixed number of DeriveDT(dt) calls covering T ; monitors every call.
// (For RK45 the call step dt only caps the substep : give it a large dt.)
result Run(const config & c, event::integrator id, Double_t dt, Double_t tol = 0.){
  event* ev = Make(c);
  ev->SetIntegrator(id);
  if(tol>0.) ev->SetTolerance(tol, tol*1e-2);
  result r;
  Double_t e0 = ev->GetNetE(), ps, ls;
  Scales(ev, ps, ls);
  TVector3 p0 = ev->GetNetMomentum(), l0 = ev->GetNetAngularMomentum();
  Long_t n = (Long_t)TMath::Nint(c.T/dt);
  for(Long_t s=0; s<n; s++){
    ev->DeriveDT(dt);
    Double_t e = TMath::Abs(ev->GetNetE()-e0)/TMath::Abs(e0);
    if(e > r.emax) r.emax = e;
    Double_t dmin = ev->GetMinDistance();
    if(dmin < r.rmin) r.rmin = dmin;
  }
  r.dp = (ev->GetNetMomentum()-p0).Mag()/ps;
  r.dl = (ev->GetNetAngularMomentum()-l0).Mag()/ls;
  for(Int_t i=0; i<ev->getNParticle(); i++) r.x.push_back(ev->getParticle(i)->GetX());
  r.nsub = (id==event::kRK45) ? ev->GetRK45NSubsteps() : n;
  delete ev;
  return r;
}

Double_t PosErr(const result & a, const result & ref){
  Double_t m = 0.;
  for(size_t i=0; i<a.x.size(); i++) m = TMath::Max(m, (a.x[i]-ref.x[i]).Mag());
  return m;
}

// ======================================================================
void Conservation(){
  printf("\n== A. conservation and convergence over a run (no closed form)\n");
  printf("   dE : max over the run |E-E0|/|E0| ; dP, dL : final drift / sum |p|, sum |x x p| ;\n");
  printf("   dx : final max position error vs the reference RK4 dt = 0.25 (fm) ; p : log2(dx(2dt)/dx(dt))\n");
  printf("   nstep : steps (RK45 : accepted substeps, the call step is %g fm/c)\n", 200.);
  struct integ{ event::integrator id; const char* name; Double_t dts[3]; Double_t order; };
  const integ ig[] = {
    {event::kTaylor,          "TAYLOR", {4., 2., 1.}, 1.},
    {event::kSymplecticEuler, "SEULER", {4., 2., 1.}, 1.},
    {event::kVelocityVerlet,  "VERLET", {4., 2., 1.}, 2.},
    {event::kRK4,             "RK4",    {16., 8., 4.}, 4.}};
  const Int_t ni = 4;
  std::vector<config> cs = Configs();

  // worst over the configurations, at the finest step, for the checks
  Double_t wE[5] = {0}, wP[5] = {0}, wL[5] = {0}, wX[5] = {0};
  Double_t pmin[4], pmax[4];
  for(Int_t k=0; k<ni; k++){ pmin[k] = 1e300; pmax[k] = -1e300; }
  Double_t refErr = 0.;

  for(size_t ic=0; ic<cs.size(); ic++){
    const config & c = cs[ic];
    result ref = Run(c, event::kRK4, 0.25);
    result chk = Run(c, event::kRK45, 1., 1e-13);   // independent estimate of the reference error
    refErr = TMath::Max(refErr, PosErr(chk, ref));
    printf("\n  %s : %d bodies, T = %.0f fm/c, smallest distance %.2f fm, |RK45(tol 1e-13, dt 1) - ref| = %.2e fm\n",
           c.name, (Int_t)c.b.size(), c.T, ref.rmin, PosErr(chk, ref));
    printf("  %-7s %6s %10s %10s %10s %10s %6s %8s\n", "integ", "dt", "dE", "dP", "dL", "dx (fm)", "p", "nstep");
    for(Int_t k=0; k<ni; k++){
      Double_t dx[3];
      for(Int_t j=0; j<3; j++){
        result r = Run(c, ig[k].id, ig[k].dts[j]);
        dx[j] = PosErr(r, ref);
        Double_t p = (j>0) ? TMath::Log(dx[j-1]/dx[j])/TMath::Log(2.) : 0.;
        printf("  %-7s %6.1f %10.2e %10s %10.2e %10.2e %6s %8ld\n", ig[k].name, ig[k].dts[j], r.emax,
               c.pCons ? Form("%10.2e", r.dp) : "-", r.dl, dx[j], (j>0) ? Form("%6.2f", p) : "", r.nsub);
        if(j==2){
          wE[k] = TMath::Max(wE[k], r.emax);
          if(c.pCons) wP[k] = TMath::Max(wP[k], r.dp);
          wL[k] = TMath::Max(wL[k], r.dl);
          wX[k] = TMath::Max(wX[k], dx[j]);
          pmin[k] = TMath::Min(pmin[k], p);
          pmax[k] = TMath::Max(pmax[k], p);
        }
      }
    }
    const Double_t tols[] = {1e-6, 1e-8, 1e-10};
    for(Int_t j=0; j<3; j++){
      result r = Run(c, event::kRK45, 200., tols[j]);
      printf("  %-7s %6s %10.2e %10s %10.2e %10.2e %6s %8ld\n", Form("RK45-%d", (Int_t)TMath::Nint(-TMath::Log10(tols[j]))),
             "tol", r.emax, c.pCons ? Form("%10.2e", r.dp) : "-", r.dl, PosErr(r, ref), "", r.nsub);
      if(j==1){
        wE[4] = TMath::Max(wE[4], r.emax);
        if(c.pCons) wP[4] = TMath::Max(wP[4], r.dp);
        wL[4] = TMath::Max(wL[4], r.dl);
        wX[4] = TMath::Max(wX[4], PosErr(r, ref));
      }
    }
  }

  printf("\n");
  Check("reference error estimate |RK45(1e-13) - RK4(dt=0.25)| (fm)", refErr, 1e-9);
  // tolerances (finest step, worst configuration) about 3-5x the measured values
  const char* names[] = {"TAYLOR dt=1", "SEULER dt=1", "VERLET dt=1", "RK4 dt=4", "RK45 tol=1e-8"};
  // measured (worst of C1-C3) : TAYLOR 2.8e-4 / 1.8e-4 / 0.21 fm ; SEULER 1.1e-4 / 6e-14 / 2.0e-3 ;
  // VERLET 4.1e-8 / 8e-14 / 2.6e-5 ; RK4 3.8e-12 / 3e-14 / 1.6e-9 ; RK45 1.6e-9 / 1.9e-10 / 3.2e-6
  // (energy / angular momentum / position). SEULER and VERLET conserve L to round-off
  // (symplectic + central pair forces) ; the 5e-13 allows round-off growth over 4e4 steps.
  const Double_t tE[] = {1e-3,  5e-4,  2e-7,  2e-11, 1e-8};
  const Double_t tX[] = {1.,    1e-2,  1e-4,  1e-8,  2e-5};
  const Double_t tL[] = {1e-3,  5e-13, 5e-13, 2e-13, 1e-9};
  for(Int_t k=0; k<5; k++){
    Check(Form("%-13s max energy error over the run", names[k]), wE[k], tE[k]);
    Check(Form("%-13s momentum drift (C1, C3)", names[k]), wP[k], 1e-13);  // pair forces are equal and opposite
    Check(Form("%-13s angular momentum drift", names[k]), wL[k], tL[k]);
    Check(Form("%-13s final position error (fm)", names[k]), wX[k], tX[k]);
  }
  for(Int_t k=0; k<ni; k++){
    CheckRange(Form("%-6s order, smallest over the configurations", ig[k].name), pmin[k], ig[k].order-0.3, ig[k].order+0.3);
    CheckRange(Form("%-6s order, largest over the configurations", ig[k].name), pmax[k], ig[k].order-0.3, ig[k].order+0.3);
  }
}

// ======================================================================
// alpha + fixed Au (+ proton at (0, D, 0)), RK45 tol 1e-12, n calls of dt
event* Decoupled(Double_t D, Bool_t withThird){
  event* ev = new event();
  ev->AddForce(&gCoulomb);
  ev->AddParticle(new EMparticle(4., 2., TVector3(-1000., 50., 0.), TVector3(0.05, 0., 0.), false, false));
  ev->AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  if(withThird) ev->AddParticle(new EMparticle(1.00728, 1., TVector3(0., D, 0.), TVector3(), false, false));
  ev->SetIntegrator(event::kRK45);
  ev->SetTolerance(1e-12, 1e-14);
  return ev;
}

void Decoupling(){
  printf("\n== B. a far third body : alpha (b = 50 fm, 0.05c, x0 = -1000) on fixed Au + proton at (0, D, 0)\n");
  printf("   theta : asymptote of the alpha-Au pair (initial) -> (final, proton ignored) ;\n");
  printf("   theta_R : exact two-body Rutherford at the asymptotic b, v (by hand)\n");
  printf("   est : first-order bound for a uniform field F = k q_a q_p / D^2 (along y) :\n");
  printf("         est = F/(m v0) [ T + |dtheta/db| Int|x(t)|dt ] -- the direct impulse F T, plus\n");
  printf("         the torque x F about the target that changes L, hence b, by Int|x|F dt / (m v0)\n");
  const Double_t kap = gCoulomb.GetK()*2.*79.;
  const Double_t dt = 10.;
  // two-body run : number of calls to the stop (receding, 1000 fm), lever-arm integral
  event* ev2 = Decoupled(0., kFALSE);
  asymptote asyIn(ev2->getParticle(0), ev2->getParticle(1), kap);
  Long_t n = 0;
  Double_t lever = 0.;   // Int |x(t)| dt  (fm^2/c)
  for(; n<1000000; n++){
    TVector3 dx = ev2->getParticle(0)->GetX(), dv = ev2->getParticle(0)->GetV();
    if(dx.Dot(dv)>0. && dx.Mag()>=1000.) break;
    Double_t x1 = TMath::Abs(dx.X());
    ev2->DeriveDT(dt);
    lever += 0.5*(x1 + TMath::Abs(ev2->getParticle(0)->GetX().X()))*dt;
  }
  const Double_t T = n*dt;
  asymptote asy2(ev2->getParticle(0), ev2->getParticle(1), kap);
  TVector3 v2 = ev2->getParticle(0)->GetV();
  // exact Rutherford, by hand
  TVector3 r0(-1000., 50., 0.), w0(0.05, 0., 0.);
  Double_t E = 0.5*4.*w0.Mag2() + kap/r0.Mag();
  Double_t vinf = TMath::Sqrt(2.*E/4.), binf = r0.Cross(w0).Mag()/vinf;
  Double_t thR = 2.*TMath::ATan(kap/(4.*vinf*vinf*binf));
  Double_t th2 = asymptote::Angle(asyIn.GetVInRel(), asy2.GetVOutRel());
  Double_t aR = kap/(4.*vinf*vinf);
  Double_t dthdb = 2.*aR/(aR*aR + binf*binf);
  printf("  two-body : T = %.0f fm/c, theta = %.9f deg, theta - theta_R = %.2e rad\n", T, th2*TMath::RadToDeg(), th2-thR);
  printf("             Int|x|dt = %.4g fm^2/c, |dtheta/db| = %.4g rad/fm : lever term / T = %.2f\n",
         lever, dthdb, dthdb*lever/T);
  delete ev2;

  const Double_t Ds[] = {1e4, 3e4, 1e5, 3e5, 1e6};
  const Int_t nD = 5;
  Double_t dev[5], est[5];
  printf("  %10s %12s %12s %12s %12s %8s\n", "D (fm)", "est (rad)", "|v3 - v2|dir", "theta-th_R", "|..|/est", "slope");
  Double_t worstRatio = 0., worstR = 0.;
  for(Int_t i=0; i<nD; i++){
    event* ev = Decoupled(Ds[i], kTRUE);
    for(Long_t s=0; s<n; s++) ev->DeriveDT(dt);
    asymptote asy3(ev->getParticle(0), ev->getParticle(1), kap);
    TVector3 v3 = ev->getParticle(0)->GetV();
    dev[i] = asymptote::Angle(v3, v2);
    est[i] = gCoulomb.GetK()*2.*1./(Ds[i]*Ds[i]*4.*0.05) * (T + dthdb*lever);
    Double_t th3 = asymptote::Angle(asyIn.GetVInRel(), asy3.GetVOutRel());
    Double_t slope = (i>0) ? TMath::Log(dev[i-1]/dev[i])/TMath::Log(Ds[i]/Ds[i-1]) : 0.;
    printf("  %10.0e %12.3e %12.3e %12.3e %12.3f %8s\n", Ds[i], est[i], dev[i], th3-thR,
           TMath::Abs(th3-thR)/est[i], (i>0) ? Form("%8.3f", slope) : "");
    worstRatio = TMath::Max(worstRatio, dev[i]/est[i]);
    worstR     = TMath::Max(worstR, TMath::Abs(th3-thR)/est[i]);
    delete ev;
  }
  Check("direction change by the third body / estimate (bound 1)", worstRatio, 1.);
  Check("|theta_3body - theta_Rutherford| / estimate (bound 1)", worstR, 1.);
  // the estimate is not a loose bound : the effect is really of that size
  Double_t least = 1e300;
  for(Int_t i=0; i<nD; i++) least = TMath::Min(least, dev[i]/est[i]);
  CheckRange("direction change / estimate, smallest (estimate is tight)", least, 0.1, 1.);
  // 1/D^2 scaling, between the first and the last D
  Double_t slope = TMath::Log(dev[0]/dev[nD-1])/TMath::Log(Ds[nD-1]/Ds[0]);
  CheckRange("deviation ~ D^-k : k from D = 1e4 to 1e6", slope, 1.9, 2.1);
}

// ======================================================================
// RK45 (tol 1e-13) from the initial state ; whenever all pairs recede and the
// smallest pair distance passes R = 1e3, 1e4, ... 1e12 fm, store the unit
// velocity of every moving body (and, for two free bodies, of the relative velocity)
struct snapR{ Double_t R; std::vector<TVector3> n; std::vector<Double_t> s; };

std::vector<snapR> Asymptotics(const config & c, Bool_t relative){
  std::vector<snapR> out;
  event* ev = Make(c);
  ev->SetIntegrator(event::kRK45);
  ev->SetTolerance(1e-13, 1e-15);
  const Int_t np = ev->getNParticle();
  Double_t Rnext = 1e3;
  //! [far_loop]
  for(Long_t g=0; g<2000000 && Rnext<=1e12*1.0001; g++){
    Double_t rmin = 1e300, vmax = 0.;
    Bool_t receding = kTRUE;
    for(Int_t i=0; i<np; i++) for(Int_t j=i+1; j<np; j++){
      particle* a = ev->getParticle(i);
      particle* b = ev->getParticle(j);
      TVector3 dx = a->GetX()-b->GetX(), dv = a->GetV()-b->GetV();
      rmin = TMath::Min(rmin, dx.Mag());
      vmax = TMath::Max(vmax, dv.Mag());
      if(dx.Dot(dv) <= 0.) receding = kFALSE;
    }
    if(receding && rmin >= Rnext){
      snapR s; s.R = rmin;
      if(relative){
        TVector3 vr = ev->getParticle(0)->GetV() - ev->getParticle(1)->GetV();
        s.n.push_back(vr.Unit()); s.s.push_back(vr.Mag());
      }else{
        for(Int_t i=0; i<np; i++){
          particle* p = ev->getParticle(i);
          if(p->IsInvincible()) continue;
          s.n.push_back(p->GetV().Unit()); s.s.push_back(p->GetV().Mag());
        }
      }
      out.push_back(s);
      Rnext *= 10.;
      continue;
    }
    ev->DeriveDT(TMath::Max(2., 0.01*rmin/vmax));
  }
  //! [far_loop]
  delete ev;
  return out;
}

void AsymptoticConvergence(){
  printf("\n== C. outgoing directions vs stopping distance R (RK45 tol 1e-13, reference R = 1e12 fm)\n");
  printf("   delta(R) = max over bodies of the angle to the R = 1e12 direction ; ds = max |speed - ref| / ref\n");
  std::vector<config> cs = Configs();
  config f2; f2.name = "alpha + fixed Au (2-body)"; f2.T = 0.; f2.pCons = kFALSE;
  f2.b.push_back({4.,   2.,  TVector3(-1000., 50., 0.), TVector3(0.05, 0., 0.), kFALSE});
  f2.b.push_back({197., 79., TVector3(), TVector3(), kTRUE});
  config fr = f2; fr.name = "alpha + free Au (2-body)";
  fr.b[1].fixed = kFALSE;

  struct acase{ const config* c; Bool_t rel; const char* label; Double_t k; };
  const acase ac[] = {
    {&f2,    kFALSE, "alpha + fixed Au, lab",        2.},
    {&fr,    kTRUE,  "alpha + free Au, relative v",  2.},
    {&fr,    kFALSE, "alpha + free Au, lab",         1.},
    {&cs[0], kFALSE, "C1 3 free, lab",               1.},
    {&cs[1], kFALSE, "C2 fixed Au + 2 alpha, lab",   1.},
    {&cs[2], kFALSE, "C3 4 free, lab",               1.}};
  const Int_t na = sizeof(ac)/sizeof(ac[0]);

  printf("  %-30s", "R (fm) :");
  for(Int_t e=3; e<=11; e++) printf(" %8.0e", TMath::Power(10., e));
  printf(" %7s %9s %10s\n", "k", "C", "R(1e-6)");
  for(Int_t a=0; a<na; a++){
    std::vector<snapR> s = Asymptotics(*ac[a].c, ac[a].rel);
    if(s.size() < 10){
      printf("  %-30s only %d stopping distances reached\n", ac[a].label, (Int_t)s.size());
      CheckRange(Form("%s : reached R = 1e12", ac[a].label), s.size(), 10, 10);
      continue;
    }
    const snapR & ref = s.back();
    std::vector<Double_t> lr, ld, ds;
    printf("  %-30s", Form("%s  delta", ac[a].label));
    for(size_t i=0; i+1<s.size(); i++){
      Double_t d = 0., dsp = 0.;
      for(size_t j=0; j<ref.n.size(); j++){
        d   = TMath::Max(d, asymptote::Angle(s[i].n[j], ref.n[j]));
        dsp = TMath::Max(dsp, TMath::Abs(s[i].s[j]-ref.s[j])/ref.s[j]);
      }
      printf(" %8.1e", d);
      lr.push_back(TMath::Log10(s[i].R)); ld.push_back(TMath::Log10(d)); ds.push_back(dsp);
    }
    // least-squares slope over R = 1e4 ... 1e8 (indices 1..5) : above the round-off
    // floor (~1e-15 rad), which a 1/R^2 case reaches near 1e9 fm
    Double_t sx=0., sy=0., sxx=0., sxy=0.; Int_t m=0;
    for(Int_t i=1; i<=5; i++){ sx+=lr[i]; sy+=ld[i]; sxx+=lr[i]*lr[i]; sxy+=lr[i]*ld[i]; m++; }
    Double_t slope = (m*sxy - sx*sy)/(m*sxx - sx*sx);
    Double_t inter = (sy - slope*sx)/m;
    Double_t k = -slope;
    Double_t C = TMath::Power(10., inter);           // delta = C / R^k  (C in rad fm^k)
    Double_t R6 = TMath::Power(C/1e-6, 1./k);        // R where delta = 1e-6 rad
    printf(" %7.3f %9.2e %10.2e\n", k, C, R6);
    printf("  %-30s", "                    ds");
    for(size_t i=0; i<ds.size(); i++) printf(" %8.1e", ds[i]);
    printf("\n");
    CheckRange(Form("%-30s rate k (1e4..1e8 fm)", ac[a].label), k, ac[a].k-0.1, ac[a].k+0.1);
  }
}

int main(){
  std::clock_t c0 = std::clock();
  printf("validation : N-body Coulomb systems (numerical evidence, no closed form)\n");
  Conservation();
  Decoupling();
  AsymptoticConvergence();
  printf("\n%d / %d checks passed  (%.1f s)\n", gCheck-gFail, gCheck, Double_t(std::clock()-c0)/CLOCKS_PER_SEC);
  return gFail==0 ? 0 : 1;
}
