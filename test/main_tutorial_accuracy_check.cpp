////////////////////////////////////////////
//
// RooParticle tutorial, part 6 : choosing the accuracy
//
// 1. Two bodies (alpha on fixed Au, b = 50 fm): every integrator, a few dt,
//    run until receding at 1000 fm. The exact answer is the asymptote of the
//    initial state; the cost is the number of force evaluations.
// 2. Three free bodies (no exact answer): the same quantity with RK45 at
//    tighter and tighter tolerance converges; the differences are the error.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"
#include "asymptote.hpp"

#include "TMath.h"
#include "TVector3.h"
#include <cstdio>

static Int_t gFail = 0;
void Expect(const char* what, Bool_t ok){
  if(!ok){ std::printf("  FAIL : %s\n", what); gFail++; }
}

//! [counting]
// Coulomb force that counts how often it is evaluated: the cost of a run.
class countingCoulomb : public coulombForce{
public:
  TVector3 Evaluate(particle* p1, particle* p2){ ncall++; return coulombForce::Evaluate(p1, p2); }
  Long_t ncall = 0;
};
//! [counting]

struct result{ Long_t steps, evals; Double_t dE, dTheta; };

result RunTwoBody(const char* integ, Double_t dt, Double_t rtol = 0.){
  countingCoulomb f;
  event ev;
  ev.AddForce(&f);
  particle* t = ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  particle* a = ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 50., 0.), TVector3(0.05, 0., 0.), false, false));
  ev.SetIntegrator(integ);
  if(rtol > 0.) ev.SetTolerance(rtol, rtol * 1e-2);
  ev.AddInspector(new inspectorP(t, a, "VRG", 0.));
  ev.AddInspector(new inspectorP(t, a, "DSTG", 1000.));

  const Double_t kappa = f.GetK() * 2. * 79.;
  asymptote exact(a, t, kappa);            // exact for this two-body problem
  Double_t E0 = ev.GetNetE();
  f.ncall = 0;
  ev.DeriveInspect(1, dt);
  asymptote got(a, t, kappa);

  result r;
  r.steps  = ev.GetDeriveN();
  r.evals  = f.ncall;
  r.dE     = TMath::Abs((ev.GetNetE() - E0) / E0);
  r.dTheta = TMath::Abs(got.GetThetaCM() - exact.GetThetaCM()) * TMath::RadToDeg();
  return r;
}

void TwoBody(){
  std::printf("== 1. alpha on fixed Au, b = 50 fm, until receding at 1000 fm ==\n");
  //! [table]
  struct cfg{ const char* integ; Double_t dt, rtol; };
  const cfg cs[] = {
    {"TAYLOR", 1., 0.}, {"TAYLOR", 0.1, 0.},
    {"SEULER", 1., 0.}, {"SEULER", 0.1, 0.},
    {"VERLET", 4., 0.}, {"VERLET", 1., 0.},
    {"RK4",    8., 0.}, {"RK4",    2., 0.},
    {"RK45", 100., 1e-6}, {"RK45", 100., 1e-10},       // dt limits the substep
    {"RK45", 2000., 1e-6}, {"RK45", 2000., 1e-8}, {"RK45", 2000., 1e-10}};
  const Int_t n = sizeof(cs) / sizeof(cs[0]);
  std::printf("%-7s %6s %8s %9s %10s %10s %15s\n", "integ", "dt", "rel tol", "steps", "force eval", "|dE/E|", "angle err (deg)");
  result r[n];
  for(Int_t i = 0; i < n; i++){
    r[i] = RunTwoBody(cs[i].integ, cs[i].dt, cs[i].rtol);
    std::printf("%-7s %6.1f %8.0e %9ld %10ld %10.1e %15.1e\n", cs[i].integ, cs[i].dt, cs[i].rtol,
                r[i].steps, r[i].evals, r[i].dE, r[i].dTheta);
  }
  //! [table]
  Expect("TAYLOR first order", r[0].dTheta / r[1].dTheta > 5. && r[0].dTheta / r[1].dTheta < 20.);
  Expect("VERLET second order", r[4].dTheta / r[5].dTheta > 10.);
  Expect("RK4 beats VERLET at fewer evaluations", r[6].dTheta < r[5].dTheta && r[6].evals < r[5].evals);
  Expect("RK45 dt=100 : the tolerance hardly matters", r[9].evals < 1.2 * r[8].evals);
  Expect("RK45 dt=2000 : tolerance controls the error", r[12].dTheta < r[11].dTheta && r[11].dTheta < r[10].dTheta);
  Expect("RK45 1e-10 cheap and accurate", r[12].dTheta < 1e-7 && r[12].evals < r[5].evals / 10);
}

//! [nbody]
// 12C at rest, an alpha and a proton, all free ; returns the alpha direction
// (deg) once all pairs recede and the closest pair is 1e4 fm apart.
Double_t AlphaAngle(Double_t rtol, Long_t & evals){
  countingCoulomb f;
  event ev;
  ev.AddForce(&f);
  ev.AddParticle(new EMparticle(12., 6., TVector3(0., 0., 0.), TVector3(), false, false));
  ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 20., 0.), TVector3(0.05, 0., 0.), false, false));
  ev.AddParticle(new EMparticle(1., 1., TVector3(-900., -60., 10.), TVector3(0.04, 0., 0.), false, false));
  ev.SetIntegrator("RK45");
  ev.SetTolerance(rtol, rtol * 1e-2);
  for(Int_t i = 0; i < 3; i++) for(Int_t j = i + 1; j < 3; j++)
    ev.AddInspector(new inspectorP(ev.getParticle(i), ev.getParticle(j), "VRG", 0.));
  ev.AddInspector(new inspectorE(&ev, "DSG", 1e4));
  ev.DeriveInspect(1, 2000.);   // large dt : RK45 chooses its own substeps inside
  evals = f.ncall;
  return TVector3(1., 0., 0.).Angle(ev.getParticle(1)->GetV()) * TMath::RadToDeg();
}
//! [nbody]

void ThreeBody(){
  std::printf("\n== 2. three free bodies, RK45, alpha direction at 1e4 fm ==\n");
  //! [nbodytable]
  const Double_t tols[] = {1e-5, 1e-7, 1e-9, 1e-11, 1e-13};
  Double_t prev = 0., diff[5];
  std::printf("%8s %12s %16s %18s\n", "rel tol", "force eval", "alpha angle (deg)", "change (deg)");
  for(Int_t i = 0; i < 5; i++){
    Long_t ev;
    Double_t a = AlphaAngle(tols[i], ev);
    diff[i] = (i > 0) ? TMath::Abs(a - prev) : 0.;
    std::printf("%8.0e %12ld %16.10f %18.1e\n", tols[i], ev, a, diff[i]);
    prev = a;
  }
  //! [nbodytable]
  Expect("three-body angle converges", diff[4] < diff[2] && diff[4] < 1e-8);
}

int main(){
  TwoBody();
  ThreeBody();
  std::printf("%s\n", gFail == 0 ? "TUTORIAL 6 OK" : "TUTORIAL 6 FAILED");
  return gFail == 0 ? 0 : 1;
}
