////////////////////////////////////////////
//
// RooParticle tutorial, part 4 : three and more bodies
//
// A. Fixed target plus two alpha particles that also repel each other.
//    Terminating with event inspectors and one pair inspector per pair;
//    energy and angular momentum about the fixed centre as the correctness
//    evidence; the effect of the alpha-alpha force.
// B. Three free bodies (12C at rest, an alpha and a proton coming in).
//    Energy, momentum and angular momentum; asymptotic velocities by
//    integrating with RK45 to larger and larger distance and watching the
//    result converge.
// C. oddevent against event once more than one body moves.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "event.hpp"
#include "oddevent.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"

#include "TMath.h"
#include "TVector3.h"
#include <cstdio>
#include <vector>

static coulombForce gCoulomb;
static Int_t gFail = 0;
void Expect(const char* what, Bool_t ok){
  if(!ok){ std::printf("  FAIL : %s\n", what); gFail++; }
}

//! [pairs]
// One "receding" inspector (radial velocity > 0) for every pair in which at
// least one particle moves. In AND mode the event goes on until all of them
// are satisfied at the same inspection. n particles give up to n(n-1)/2 pairs.
void AddRecedingInspectors(event & ev){
  for(Int_t i = 0; i < ev.getNParticle(); i++){
    for(Int_t j = i + 1; j < ev.getNParticle(); j++){
      particle* pi = ev.getParticle(i);
      particle* pj = ev.getParticle(j);
      if(pi->IsInvincible() && pj->IsInvincible()) continue;
      ev.AddInspector(new inspectorP(pi, pj, "VRG", 0.));
    }
  }
}
//! [pairs]

//! [termination]
// Terminating conditions that make sense for any number of bodies.
void AddManyBodyStop(event & ev, Double_t pkl){
  AddRecedingInspectors(ev);                        // every pair flies apart
  ev.AddInspector(new inspectorE(&ev, "PKL", pkl)); // |PE / KE| < pkl : interaction small
  inspector* drift = new inspectorE(&ev, "EDG", 1e-6);  // numerical failure
  drift->SetAbort();
  ev.AddInspector(drift);
  ev.SetMaxInspectSteps(2000000);
}
//! [termination]

// target and alpha 1 (and alpha 2 if both)
void BuildA(event & ev, Bool_t a1, Bool_t a2){
  ev.AddForce(&gCoulomb);
  ev.SetIntegrator("VERLET");
  ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  if(a1) ev.AddParticle(new EMparticle(4., 2., TVector3(-1000.,  40.,  0.), TVector3(0.05, 0., 0.), false, false));
  if(a2) ev.AddParticle(new EMparticle(4., 2., TVector3(-1060., -25., 15.), TVector3(0.05, 0., 0.), false, false));
}

Double_t Deg(TVector3 a, TVector3 b){ return a.Angle(b) * TMath::RadToDeg(); }

void PartA(){
  std::printf("== A. fixed Au + two alphas (VERLET, dt = 1 fm/c) ==\n");
  //! [fixed]
  event ev;
  BuildA(ev, kTRUE, kTRUE);
  AddManyBodyStop(ev, 0.01);
  std::printf("inspectors          : %d (3 VRG pairs, PKL, abort EDG)\n", ev.getNInspector());

  const TVector3 origin = ev.getParticle(0)->GetX();   // the fixed centre
  const Double_t E0 = ev.GetNetE();
  const TVector3 L0 = ev.GetNetAngularMomentum(origin);

  ev.DeriveInspect();
  std::printf("stop                : %s after %ld steps\n", ev.GetStopReasonName().Data(), (long) ev.GetDeriveN());
  std::printf("|dE/E|              : %.2e\n", TMath::Abs((ev.GetNetE() - E0) / E0));
  std::printf("|dL|/|L| (about Au) : %.2e\n", (ev.GetNetAngularMomentum(origin) - L0).Mag() / L0.Mag());
  std::printf("|P| (not conserved, the fixed centre takes the recoil) : %.4f -> %.4f amu c\n",
              8. * 0.05, ev.GetNetMomentum().Mag());
  //! [fixed]
  Expect("A ends by inspectors", ev.GetStopReason() == event::kStopInspector);
  Expect("A energy", TMath::Abs((ev.GetNetE() - E0) / E0) < 1e-8);
  Expect("A angular momentum", (ev.GetNetAngularMomentum(origin) - L0).Mag() / L0.Mag() < 1e-10);

  //! [alone]
  // The same alphas, each one alone with the target: the difference is the
  // effect of the alpha-alpha force.
  for(Int_t k = 1; k <= 2; k++){
    event one;
    BuildA(one, k == 1, k == 2);
    AddManyBodyStop(one, 0.01);
    one.DeriveInspect();
    TVector3 v3 = ev.getParticle(k)->GetV(), v2 = one.getParticle(1)->GetV();
    std::printf("alpha %d : angle %.4f deg with the other alpha, %.4f deg alone, |v| %.6f vs %.6f c\n",
                k, Deg(TVector3(1., 0., 0.), v3), Deg(TVector3(1., 0., 0.), v2), v3.Mag(), v2.Mag());
    Expect("A alpha-alpha force visible", Deg(v3, v2) > 0.1);
  }
  //! [alone]
}

// 12C at rest, alpha and proton coming in ; all free
void BuildB(event & ev){
  ev.AddForce(&gCoulomb);
  ev.AddParticle(new EMparticle(12., 6., TVector3(0., 0., 0.), TVector3(), false, false));
  ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 20., 0.), TVector3(0.05, 0., 0.), false, false));
  ev.AddParticle(new EMparticle(1., 1., TVector3(-900., -60., 10.), TVector3(0.04, 0., 0.), false, false));
}

void PartB(){
  std::printf("\n== B. three free bodies: 12C at rest, alpha, proton (RK45) ==\n");
  //! [free]
  event ev;
  BuildB(ev);
  ev.SetIntegrator("RK45");
  ev.SetTolerance(1e-10, 1e-12);

  const Double_t E0 = ev.GetNetE();
  const TVector3 P0 = ev.GetNetMomentum();
  const TVector3 L0 = ev.GetNetAngularMomentum();
  Double_t pscale = 0.;                       // scale for the momentum error
  for(Int_t i = 0; i < ev.getNParticle(); i++) pscale += ev.getParticle(i)->GetMomentum().Mag();

  AddRecedingInspectors(ev);
  inspector* far = new inspectorE(&ev, "DSG", 1e3);   // smallest pair distance > R
  ev.AddInspector(far);
  //! [free]

  //! [converge]
  // Integrate to R = 1e3, 1e4, ... 1e7 fm (the same event continues) and
  // print the velocities there. The change from one R to the next shrinks
  // like 1/R: that is the convergence test. |PE/E| says how much energy is
  // still in the interaction.
  std::printf("%9s %12s %12s %12s %10s %10s %10s\n", "R (fm)", "|v C| (c)", "|v alpha|", "|v p|",
              "angle a", "angle p", "|PE/E|");
  std::vector<TVector3> prev;
  Double_t lastChange = 0.;
  for(Double_t R = 1e3; R <= 1e7 * 1.0001; R *= 10.){
    far->SetValue(R);
    ev.DeriveInspect(1, R / 100.);          // dt grows with R; RK45 subdivides as needed
    std::vector<TVector3> v;
    for(Int_t i = 0; i < 3; i++) v.push_back(ev.getParticle(i)->GetV());
    Double_t change = 0.;
    if(!prev.empty()) for(Int_t i = 0; i < 3; i++) change = TMath::Max(change, (v[i] - prev[i]).Mag() / v[i].Mag());
    std::printf("%9.0e %12.9f %12.9f %12.9f %10.6f %10.6f %10.2e   max rel. change %.1e\n",
                R, v[0].Mag(), v[1].Mag(), v[2].Mag(),
                Deg(TVector3(1., 0., 0.), v[1]), Deg(TVector3(1., 0., 0.), v[2]),
                TMath::Abs(ev.GetNetPE() / ev.GetNetE()), change);
    if(R >= 1e7 * 0.9999){
      // the tail is ~ 1/R : v_inf ~ v(R) + (v(R) - v(R/10)) / 9   (Richardson)
      for(Int_t i = 0; i < 3; i++){
        TVector3 vinf = v[i] + (v[i] - prev[i]) * (1. / 9.);
        std::printf("extrapolated |v_inf| of particle %d : %.9f c (correction %.1e)\n",
                    i, vinf.Mag(), (vinf - v[i]).Mag() / v[i].Mag());
      }
    }
    prev = v;
    lastChange = change;
  }
  //! [converge]

  //! [conserved]
  std::printf("|dE/E|              : %.2e\n", TMath::Abs((ev.GetNetE() - E0) / E0));
  std::printf("|dP| / sum |p|      : %.2e\n", (ev.GetNetMomentum() - P0).Mag() / pscale);
  std::printf("|dL| / |L|          : %.2e\n", (ev.GetNetAngularMomentum() - L0).Mag() / L0.Mag());
  std::printf("RK45                : %ld steps, %ld substeps, %ld rejected\n",
              (long) ev.GetDeriveN(), (long) ev.GetRK45NSubsteps(), (long) ev.GetRK45NRejected());
  //! [conserved]
  Expect("B ends by inspectors", ev.GetStopReason() == event::kStopInspector);
  Expect("B energy", TMath::Abs((ev.GetNetE() - E0) / E0) < 1e-8);
  Expect("B momentum", (ev.GetNetMomentum() - P0).Mag() / pscale < 1e-12);
  Expect("B angular momentum", (ev.GetNetAngularMomentum() - L0).Mag() / L0.Mag() < 1e-8);
  Expect("B converged to 2e-5", lastChange < 2e-5);
}

void PartC(){
  std::printf("\n== C. oddevent vs event (TAYLOR), dt = 1 fm/c, 60000 steps ==\n");
  //! [odd]
  // One moving particle on a fixed target: identical, bit for bit.
  event  e1;  BuildA(e1, kTRUE, kFALSE);  e1.SetIntegrator("TAYLOR");
  oddevent o1; BuildA(o1, kTRUE, kFALSE);
  e1.DeriveDTN(60000, 1.);
  o1.DeriveDTN(60000, 1.);
  std::printf("one moving body     : |x_event - x_oddevent| = %.1e fm\n",
              (e1.getParticle(1)->GetX() - o1.getParticle(1)->GetX()).Mag());

  // Three free bodies: oddevent moves each particle right after its own force
  // is computed, so action and reaction are taken at different states.
  event  e3;  BuildB(e3);  e3.SetIntegrator("TAYLOR");
  oddevent o3; BuildB(o3);
  const TVector3 P0 = e3.GetNetMomentum();
  e3.DeriveDTN(60000, 1.);
  o3.DeriveDTN(60000, 1.);
  std::printf("three free bodies   : |dP| event %.1e, oddevent %.1e amu c ; |x_event - x_oddevent| (alpha) = %.2e fm\n",
              (e3.GetNetMomentum() - P0).Mag(), (o3.GetNetMomentum() - P0).Mag(),
              (e3.getParticle(1)->GetX() - o3.getParticle(1)->GetX()).Mag());
  //! [odd]
  Expect("C identical for one moving body", e1.getParticle(1)->GetX() == o1.getParticle(1)->GetX());
  Expect("C event conserves P", (e3.GetNetMomentum() - P0).Mag() < 1e-14);
  Expect("C oddevent does not", (o3.GetNetMomentum() - P0).Mag() > 1e-9);
}

int main(){
  PartA();
  PartB();
  PartC();
  std::printf("%s\n", gFail == 0 ? "TUTORIAL 4 OK" : "TUTORIAL 4 FAILED");
  return gFail == 0 ? 0 : 1;
}
