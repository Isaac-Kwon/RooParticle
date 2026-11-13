////////////////////////////////////////////
//
// RooParticle tutorial, part 7 : beyond a fixed Coulomb centre
//
// 1. Recoiling (free) target.
// 2. A short-range nuclear force added to Coulomb (Woods-Saxon, Yukawa).
// 3. Softening of an attractive Coulomb force.
// 4. Relativistic mode.
// 5. Writing your own force (screened Coulomb).
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
#include <cmath>

static Int_t gFail = 0;
void Expect(const char* what, Bool_t ok){
  if(!ok){ std::printf("  FAIL : %s\n", what); gFail++; }
}
Double_t Deg(Double_t rad){ return rad * TMath::RadToDeg(); }

//! [myforce]
// Screened Coulomb (Bohr): U = k q1 q2 exp(-r/lambda) / r.
// Evaluate(p1, p2) is the force on p1 by p2 [amu c^2 / fm] and must be
// -grad_{x1} Potential(p1, p2), otherwise the energy monitors lie.
class screenedCoulomb : public force{
public:
  screenedCoulomb(Double_t lambda_) : force(), lambda(lambda_){}
  TVector3 Evaluate(particle* p1, particle* p2){
    TVector3 d = p1->GetX() - p2->GetX();
    Double_t r = d.Mag();
    Double_t U = k * p1->GetQ() * p2->GetQ() * std::exp(-r / lambda) / r;
    Double_t dUdr = -U * (1. / r + 1. / lambda);
    return d * (-dUdr / r);                 // F = -dU/dr * r_hat
  }
  Double_t Potential(particle* p1, particle* p2){
    Double_t r = (p1->GetX() - p2->GetX()).Mag();
    return k * p1->GetQ() * p2->GetQ() * std::exp(-r / lambda) / r;
  }
private:
  Double_t lambda;                          // screening length [fm]
  Double_t k = 0.00154629;                  // same constant as coulombForce
};
//! [myforce]

// stop when receding and far
void AddStop(event & ev, Double_t r){
  ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "VRG", 0.));
  ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", r));
}

void Recoil(){
  std::printf("== 1. alpha on 12C, fixed vs free target (b = 20 fm, v = 0.05 c) ==\n");
  //! [recoil]
  coulombForce coulomb;
  const Double_t kappa = coulomb.GetK() * 2. * 6.;
  for(Int_t fixed = 1; fixed >= 0; fixed--){
    event ev;
    ev.AddForce(&coulomb);
    ev.SetIntegrator("VERLET");
    particle* c = ev.AddParticle(new EMparticle(12., 6., TVector3(), TVector3(), fixed == 1, false));
    particle* a = ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 20., 0.), TVector3(0.05, 0., 0.), false, false));
    AddStop(ev, 1000.);
    TVector3 P0 = ev.GetNetMomentum();
    ev.DeriveInspect();
    asymptote asy(a, c, kappa);
    std::printf("%s : thetaCM %.4f deg, alpha lab %.4f deg, 12C lab %.4f deg, |dP|/|P| %.1e\n",
                fixed ? "fixed 12C" : "free 12C ", Deg(asy.GetThetaCM()), Deg(asy.GetThetaLab1()),
                Deg(asy.GetThetaLab2()), (ev.GetNetMomentum() - P0).Mag() / P0.Mag());
    if(!fixed){
      Expect("lab angle = LabAngle(thetaCM)",
             TMath::Abs(asy.GetThetaLab1() - asymptote::LabAngle(asy.GetThetaCM(), 4., 12.)) < 1e-3);
      Expect("momentum conserved", (ev.GetNetMomentum() - P0).Mag() / P0.Mag() < 1e-12);
    }
  }
  //! [recoil]
}

Double_t NuclearAngle(Double_t b, force* extra, Double_t & rmin){
  coulombForce coulomb;
  event ev;
  ev.AddForce(&coulomb);
  if(extra) ev.AddForce(extra);
  ev.SetIntegrator("RK45");
  ev.SetTolerance(1e-10, 1e-12);
  ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  particle* a = ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., b, 0.), TVector3(0.12, 0., 0.), false, false));
  AddStop(ev, 1000.);
  rmin = 1e300;
  while(!ev.Inspect() && ev.GetDeriveN() < 1000000){   // hand-made loop to watch the distance
    ev.DeriveDT(0.5);
    rmin = TMath::Min(rmin, ev.GetMinDistance());
  }
  return TVector3(1., 0., 0.).Angle(a->GetV());
}

void Nuclear(){
  std::printf("\n== 2. alpha on fixed Au at v = 0.12 c (26.8 MeV): Coulomb + nuclear ==\n");
  //! [nuclear]
  woodsSaxonForce ws;                  // V0 = 50 MeV, R = 8.886 fm, a = 0.65 fm (defaults)
  yukawaForce     yu;                  // g = 0.1 amu c^2 fm, range 1.4138 fm (illustrative)
  std::printf("%6s %14s %14s %14s %10s\n", "b (fm)", "Coulomb (deg)", "+WoodsSaxon", "+Yukawa", "rmin (fm)");
  Double_t angC[2], angW[2], angY[2], rmin;
  const Double_t bs[2] = {2., 30.};
  for(Int_t i = 0; i < 2; i++){
    angC[i] = Deg(NuclearAngle(bs[i], nullptr, rmin));
    Double_t r0 = rmin;
    angW[i] = Deg(NuclearAngle(bs[i], &ws, rmin));
    angY[i] = Deg(NuclearAngle(bs[i], &yu, rmin));
    std::printf("%6.1f %14.4f %14.4f %14.4f %10.2f\n", bs[i], angC[i], angW[i], angY[i], r0);
  }
  //! [nuclear]
  Expect("WS changes the b=2 angle", TMath::Abs(angW[0] - angC[0]) > 1.);
  Expect("WS does not change b=30", TMath::Abs(angW[1] - angC[1]) < 1e-3);
}

void Softening(){
  std::printf("\n== 3. attractive Coulomb, head-on: q = -1 on fixed q = +79 ==\n");
  //! [soft]
  const char* integ[3] = {"VERLET", "RK45", "RK45"};
  for(Int_t k = 0; k < 3; k++){
    const Bool_t soft = (k == 2);
    coulombForce coulomb;
    if(soft) coulomb.SetSoftening(2.);          // U = k q1 q2 / sqrt(r^2 + eps^2), eps = 2 fm
    event ev;
    ev.AddForce(&coulomb);
    ev.SetIntegrator(integ[k]);
    ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
    particle* p = ev.AddParticle(new EMparticle(1., -1., TVector3(-1000., 0., 0.), TVector3(0.05, 0., 0.), false, false));
    AddStop(ev, 1000.);
    inspector* crash = new inspectorE(&ev, "EDG", 1e-6);   // numerical failure
    crash->SetAbort();
    ev.AddInspector(crash);
    Double_t E0 = ev.GetNetE();
    ev.DeriveInspect(1, 1.);
    std::printf("%-6s softening %.0f fm : %-9s after %6ld steps, x = %8.1f fm, |v| = %.6f c, |dE/E| = %.1e\n",
                ev.GetIntegratorName().Data(), coulomb.GetSoftening(), ev.GetStopReasonName().Data(), (long) ev.GetDeriveN(),
                p->GetX().X(), p->GetV().Mag(), TMath::Abs((ev.GetNetE() - E0) / E0));
    if(k == 0) Expect("unsoftened VERLET aborts", ev.GetStopReason() == event::kStopAbort);
    if(k == 1) Expect("unsoftened RK45 gives up", ev.GetStopReason() == event::kStopFailed && ev.HasIntegratorFailed());
    if(k == 2) Expect("softened passes through", ev.GetStopReason() == event::kStopInspector && p->GetX().X() > 999.);
  }
  //! [soft]
}

void Relativistic(){
  std::printf("\n== 4. alpha on fixed Au at v = 0.5 c, b = 20 fm (RK4, dt = 0.25 fm/c) ==\n");
  //! [rel]
  coulombForce coulomb;
  for(Int_t rel = 0; rel <= 1; rel++){
    event ev;
    ev.AddForce(&coulomb);
    ev.SetIntegrator("RK4");             // TAYLOR (and oddevent) would ignore SetRelativistic
    ev.SetRelativistic(rel == 1);
    ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
    particle* a = ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 20., 0.), TVector3(0.5, 0., 0.), false, false));
    AddStop(ev, 1000.);
    Double_t E0 = ev.GetNetE(), E0rel = ev.GetNetERel();
    ev.DeriveInspect(1, 0.25);
    std::printf("%-17s : angle %.4f deg, |v| %.6f c, |dE/E| Newton %.1e, relativistic %.1e\n",
                rel ? "relativistic" : "non-relativistic", Deg(TVector3(1., 0., 0.).Angle(a->GetV())),
                a->GetV().Mag(), TMath::Abs((ev.GetNetE() - E0) / E0), TMath::Abs((ev.GetNetERel() - E0rel) / E0rel));
    if(rel) Expect("relativistic energy conserved", TMath::Abs((ev.GetNetERel() - E0rel) / E0rel) < 1e-9);
    else    Expect("Newtonian energy conserved", TMath::Abs((ev.GetNetE() - E0) / E0) < 1e-9);
  }
  //! [rel]
}

void MyForce(){
  std::printf("\n== 5. own force: screened Coulomb, lambda = 100 fm ==\n");
  //! [myforceuse]
  screenedCoulomb sc(100.);
  event ev;
  ev.AddForce(&sc);
  ev.SetIntegrator("VERLET");
  particle* t = ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  particle* a = ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 50., 0.), TVector3(0.05, 0., 0.), false, false));

  // check 1: force = -grad U by central differences at the start
  TVector3 x0 = a->GetX(), F = sc.Evaluate(a, t), Fn;
  const Double_t h = 1e-3;
  for(Int_t i = 0; i < 3; i++){
    TVector3 e; e(i) = h;
    a->applyX(x0 + e); Double_t up = sc.Potential(a, t);
    a->applyX(x0 - e); Double_t um = sc.Potential(a, t);
    Fn(i) = -(up - um) / (2. * h);
  }
  a->applyX(x0);
  // check 2: action = -reaction (needed by the default pair-symmetric evaluation)
  Double_t sym = (sc.Evaluate(a, t) + sc.Evaluate(t, a)).Mag() / F.Mag();

  AddStop(ev, 1000.);
  Double_t E0 = ev.GetNetE();
  ev.DeriveInspect();
  std::printf("|F - (-grad U)| / |F| = %.1e, |F12 + F21| / |F| = %.1e\n", (F - Fn).Mag() / F.Mag(), sym);
  std::printf("angle %.4f deg (unscreened: 51.01), |dE/E| = %.1e\n",
              Deg(TVector3(1., 0., 0.).Angle(a->GetV())), TMath::Abs((ev.GetNetE() - E0) / E0));
  //! [myforceuse]
  Expect("own force gradient", (F - Fn).Mag() / F.Mag() < 1e-6);
  Expect("own force symmetric", sym < 1e-14);
  Expect("own force energy", TMath::Abs((ev.GetNetE() - E0) / E0) < 1e-9);
}

int main(){
  Recoil();
  Nuclear();
  Softening();
  Relativistic();
  MyForce();
  std::printf("%s\n", gFail == 0 ? "TUTORIAL 7 OK" : "TUTORIAL 7 FAILED");
  return gFail == 0 ? 0 : 1;
}
