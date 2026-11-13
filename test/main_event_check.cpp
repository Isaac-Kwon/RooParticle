////////////////////////////////////////////
//
// RooParticle Test Application
// event check : adaptive integrator, relativistic option,
// conserved quantities, pair symmetry and ownership.
//
// 1. RK45 against a fine-step RK4 reference, for several tolerances.
// 2. Relativistic dynamics : (gamma-1) m0 + PE conserved at v ~ 0.5c,
//    same result as the non-relativistic one at v ~ 0.001c.
// 3. Total momentum (and angular momentum) of two free charged particles.
// 4. Pair-symmetric force evaluation on / off.
// 5. Events with inspectors created, copied, pre-derived and deleted in a loop.
//
// No argument. Prints a table; exit code is non-zero when a check fails.
// Rows marked "info" are printed for information only.
//
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "event.hpp"
#include "oddevent.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"

#include "iostream"
#include <cstdio>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"

static coulombForce gCoulomb;
static Int_t gFail = 0;

// Coulomb force counting its evaluations
class countingCoulomb : public coulombForce{
public:
  countingCoulomb():coulombForce(){}
  TVector3 Evaluate(particle* p1_, particle* p2_){ncall++; return coulombForce::Evaluate(p1_, p2_);}
  Long_t ncall = 0;
};

void Header(const char* title){
  printf("\n== %s ==\n", title);
  printf("%-58s %13s %11s  %s\n", "item", "value", "limit", "result");
}

// value < limit is required (kind 0), or only printed (info)
void Check(const char* name, Double_t value, Double_t limit){
  Bool_t ok = (value < limit); // false for NaN
  printf("%-58s %13.4e %11.1e  %s\n", name, value, limit, ok ? "ok" : "FAIL");
  if(!ok) gFail++;
}
void CheckTrue(const char* name, Bool_t ok, const char* what=""){
  printf("%-58s %25s  %s\n", name, what, ok ? "ok" : "FAIL");
  if(!ok) gFail++;
}
void Info(const char* name, Double_t value){
  printf("%-58s %13.4e %11s  %s\n", name, value, "", "info");
}

Double_t AngleDeg(TVector3 a, TVector3 b){
  return TMath::ATan2(a.Cross(b).Mag(), a.Dot(b)) * TMath::RadToDeg();
}

// Target (invincible) at origin and one incident particle. The event owns the particles.
void BuildFixed(event & ev, Double_t b, Double_t x0, TVector3 v,
                Double_t mT=206, Double_t qT=79, Double_t m=4, Double_t q=2){
  ev.AddForce(&gCoulomb);
  ev.AddParticle(new EMparticle(mT, qT, TVector3(0.,0.,0.), TVector3(0.,0.,0.), true, false));
  ev.AddParticle(new EMparticle(m, q, TVector3(x0, b, 0.), v, false, false));
}

// Two free charged particles
void BuildFree(event & ev, force * f = &gCoulomb){
  ev.AddForce(f);
  ev.AddParticle(new EMparticle(4,  2, TVector3(-500., 20., 0.), TVector3(0.03, 0., 0.), false, false));
  ev.AddParticle(new EMparticle(12, 6, TVector3(500., -20., 0.), TVector3(-0.02, 0.005, 0.), false, false));
}

//------------------------------------------------------------------------------
// 1. RK45
//------------------------------------------------------------------------------
void CheckRK45(){
  Header("RK45 against RK4 (alpha on Au, b = 50 fm, T = 40000 fm/c)");
  const Double_t b = 50., T = 40000., dt = 8.;
  const TVector3 vin(0.05, 0., 0.);

  event ref; BuildFixed(ref, b, -1000., vin);
  ref.SetIntegrator(event::kRK4);
  ref.DeriveDTN(TMath::Nint(T/0.125), 0.125);
  TVector3 xref = ref.getParticle(1)->GetX();
  TVector3 vref = ref.getParticle(1)->GetV();
  Info("reference RK4 dt=0.125 : scattering angle (deg)", AngleDeg(vin, vref));

  // dt = 8 : the substep is limited by dt, the tolerance hardly matters.
  // dt = 1000 : the substep is chosen by the tolerance.
  const Int_t ntol = 4;
  const Double_t tols[ntol] = {1e-6, 1e-8, 1e-10, 1e-6};
  const Double_t dtt[ntol]  = {1000., 1000., 1000., dt};
  Double_t errA[ntol], errX[ntol];
  printf("%-8s %6s %10s %10s %14s %14s %10s %8s\n",
         "rel tol", "dt", "substeps", "rejected", "angle err(deg)", "pos err (fm)", "|dE/E|", "steps");
  Bool_t timeOK = kTRUE;
  for(Int_t k=0; k<ntol; k++){
    event ev; BuildFixed(ev, b, -1000., vin);
    ev.SetIntegrator("RK45");
    ev.SetTolerance(tols[k], tols[k]*1e-2);
    Double_t e0 = ev.GetNetE();
    Int_t nstep = TMath::Nint(T/dtt[k]);
    ev.DeriveDTN(nstep, dtt[k]);
    errA[k] = AngleDeg(vref, ev.getParticle(1)->GetV());
    errX[k] = (ev.getParticle(1)->GetX() - xref).Mag();
    printf("%-8.0e %6.0f %10ld %10ld %14.3e %14.3e %10.2e %8ld\n", tols[k], dtt[k],
           (long) ev.GetRK45NSubsteps(), (long) ev.GetRK45NRejected(), errA[k], errX[k],
           TMath::Abs((ev.GetNetE()-e0)/e0), (long) ev.GetDeriveN());
    // DeriveDT(dt) : exactly dt of time and one counted step per call
    if(ev.GetDeriveN()!=nstep || ev.GetTime()!=T || ev.getParticle(1)->GetTime()!=T) timeOK = kFALSE;
  }
  CheckTrue("RK45 : N calls of DeriveDT(dt) give N steps and time N*dt", timeOK);
  Check("RK45 dt=1000 tol 1e-6  : angle error (deg)", errA[0], 1e-3);
  Check("RK45 dt=1000 tol 1e-10 : angle error (deg)", errA[2], 1e-7);
  Check("RK45 dt=1000 tol 1e-10 : position error (fm)", errX[2], 1e-5);
  CheckTrue("RK45 dt=1000 : angle error shrinks with the tolerance",
            errA[2]<errA[1] && errA[1]<errA[0]);
  CheckTrue("RK45 dt=1000 : position error shrinks with the tolerance",
            errX[2]<errX[1] && errX[1]<errX[0]);
  Check("RK45 dt=8 tol 1e-6 : angle error (deg)", errA[3], 1e-9);

  { // one call over the whole time
    event ev; BuildFixed(ev, b, -1000., vin);
    ev.SetIntegrator(event::kRK45);
    ev.SetTolerance(1e-10, 1e-12);
    ev.DeriveDT(T);
    Check("RK45 one DeriveDT(40000), tol 1e-10 : angle error (deg)",
          AngleDeg(vref, ev.getParticle(1)->GetV()), 1e-7);
    CheckTrue("RK45 one DeriveDT(40000) : 1 step, time 40000",
              ev.GetDeriveN()==1 && ev.GetTime()==T,
              TString::Format("%ld substeps", (long) ev.GetRK45NSubsteps()).Data());
  }
  { // the substep suggested by the last call is reused
    event ev; BuildFixed(ev, b, -1000., vin);
    ev.SetIntegrator(event::kRK45);
    ev.DeriveDT(dt);
    Double_t h1 = ev.GetRK45Substep();
    Long_t   n1 = ev.GetRK45NSubsteps() + ev.GetRK45NRejected();
    ev.DeriveDT(1000.);
    Long_t   n2 = ev.GetRK45NSubsteps() + ev.GetRK45NRejected() - n1;
    Info("RK45 : substep suggested after DeriveDT(8) (fm/c)", h1);
    CheckTrue("RK45 : suggestion kept above dt, then reused",
              h1>dt && n2<1000./dt, TString::Format("%ld substeps for 1000", (long) n2).Data());
  }
  { // under DeriveInspect, and copied with its settings
    event evT; BuildFixed(evT, b, -1000., vin);
    evT.SetIntegrator("RK45");
    evT.SetTolerance(1e-9, 1e-11);
    evT.AddInspector(new inspectorP(evT.getParticle(0), evT.getParticle(1), "DEG", 1));
    evT.AddInspector(new inspectorE(&evT, "TMG", 30000));
    event ev(evT);
    Bool_t rtn = ev.DeriveInspect(1, dt);
    CheckTrue("RK45 copy : DeriveInspect(1, 8) ends by the inspectors",
              rtn && ev.GetStopReason()==event::kStopInspector,
              TString::Format("time %.0f", ev.GetTime()).Data());
    CheckTrue("RK45 copy : integrator and tolerances copied",
              ev.GetIntegrator()==event::kRK45 && ev.GetRelTolerance()==1e-9 && ev.GetAbsTolerance()==1e-11);
  }
}

//------------------------------------------------------------------------------
// 2. Relativistic
//------------------------------------------------------------------------------
// largest |E(t) - E(0)| over the run, relative to the kinetic energy exchanged
Double_t RelEnergyError(event::integrator integ, Double_t dt, Bool_t rel, Double_t & dKE){
  const Double_t T = 4000.;
  event ev; BuildFixed(ev, 50., -1000., TVector3(0.5, 0., 0.));
  ev.SetIntegrator(integ);
  ev.SetRelativistic(rel);
  Double_t e0 = ev.GetNetERel();
  Double_t ke0 = ev.GetNetKERel(), kemin = ke0;
  Double_t emax = 0.;
  Int_t n = TMath::Nint(T/dt);
  for(Int_t i=0; i<n; i++){
    ev.DeriveDT(dt);
    emax  = TMath::Max(emax, TMath::Abs(ev.GetNetERel()-e0));
    kemin = TMath::Min(kemin, ev.GetNetKERel());
  }
  dKE = ke0 - kemin;
  return emax/dKE;
}

void CheckRelativistic(){
  Header("Relativistic (alpha on Au, b = 50 fm, v = 0.5c)");
  Double_t dKE = 0.;
  Double_t r;

  r = RelEnergyError(event::kRK4, 1., kTRUE, dKE);
  Info("kinetic energy exchanged (amu), RK4", dKE);
  Check("RK4    dt=1   : max|dE_rel| / dKE", r, 1e-6);
  r = RelEnergyError(event::kRK45, 10., kTRUE, dKE);
  Check("RK45   dt=10  : max|dE_rel| / dKE", r, 1e-5);
  r = RelEnergyError(event::kVelocityVerlet, 0.5, kTRUE, dKE);
  Check("VERLET dt=0.5 : max|dE_rel| / dKE", r, 1e-3);
  r = RelEnergyError(event::kSymplecticEuler, 0.5, kTRUE, dKE);
  Check("SEULER dt=0.5 : max|dE_rel| / dKE", r, 1e-1);
  r = RelEnergyError(event::kRK4, 1., kFALSE, dKE);
  Info("RK4 non-relativistic at 0.5c : max|dE_rel| / dKE", r);

  { // copy keeps the flag ; TAYLOR ignores it
    event evT; BuildFixed(evT, 50., -1000., TVector3(0.5, 0., 0.));
    evT.SetIntegrator("RK4");
    evT.SetRelativistic();
    evT.SetPairSymmetric(kFALSE);
    event ev(evT);
    CheckTrue("copy keeps SetRelativistic and SetPairSymmetric",
              ev.IsRelativistic() && !ev.IsPairSymmetric());

    event a; BuildFixed(a, 50., -1000., TVector3(0.5, 0., 0.));
    event c; BuildFixed(c, 50., -1000., TVector3(0.5, 0., 0.));
    a.SetRelativistic(kTRUE);
    a.DeriveDTN(2000, 1.);
    c.DeriveDTN(2000, 1.);
    CheckTrue("TAYLOR ignores SetRelativistic (identical result)",
              a.getParticle(1)->GetX()==c.getParticle(1)->GetX() &&
              a.getParticle(1)->GetV()==c.getParticle(1)->GetV());
  }

  Header("Relativistic, low speed (p on p fixed, b = 3000 fm, v = 0.001c)");
  {
    const Double_t T = 2e8, dt = 1e4;
    const TVector3 vin(0.001, 0., 0.);
    event er; BuildFixed(er, 3000., -1e5, vin, 1., 1., 1., 1.);
    event en; BuildFixed(en, 3000., -1e5, vin, 1., 1., 1., 1.);
    er.SetIntegrator("RK4"); er.SetRelativistic(kTRUE);
    en.SetIntegrator("RK4");
    Double_t e0r = er.GetNetERel(), e0n = en.GetNetE();
    Check("|KE_rel - KE| / KE at v = 0.001c (about 3/4 v^2)",
          TMath::Abs(er.GetNetKERel()-er.GetNetKE())/er.GetNetKE(), 1e-5);
    Int_t n = TMath::Nint(T/dt);
    er.DeriveDTN(n, dt);
    en.DeriveDTN(n, dt);
    TVector3 vr = er.getParticle(1)->GetV(), vn = en.getParticle(1)->GetV();
    Info("scattering angle (deg), non-relativistic", AngleDeg(vin, vn));
    Info("scattering angle (deg), relativistic", AngleDeg(vin, vr));
    Check("|v_rel - v_nonrel| / |v| at the end", (vr-vn).Mag()/vn.Mag(), 1e-4);
    Check("|x_rel - x_nonrel| / |x| at the end",
          (er.getParticle(1)->GetX()-en.getParticle(1)->GetX()).Mag()/en.getParticle(1)->GetX().Mag(), 1e-4);
    Check("relativistic |dE_rel / E_rel|", TMath::Abs((er.GetNetERel()-e0r)/e0r), 1e-8);
    Info("non-relativistic |dE / E|", TMath::Abs((en.GetNetE()-e0n)/e0n));
  }
}

//------------------------------------------------------------------------------
// 3. Conserved quantities
//------------------------------------------------------------------------------
void Conservation(event & ev, Long_t n, Double_t dt, Double_t & dP, Double_t & dL){
  // scale : sum of |p| of the particles, |L| at the start
  Double_t pscale = 0.;
  for(Int_t i=0; i<ev.getNParticle(); i++){
    particle* p = ev.getParticle(i);
    pscale += p->GetM0() * p->GetV().Mag();
  }
  TVector3 P0 = ev.GetNetMomentum();
  TVector3 L0 = ev.GetNetAngularMomentum();
  ev.DeriveDTN(n, dt);
  dP = (ev.GetNetMomentum() - P0).Mag() / pscale;
  dL = (ev.GetNetAngularMomentum() - L0).Mag() / L0.Mag();
}

void CheckConservation(){
  Header("Two free charged particles (m 4, q 2 and m 12, q 6), T = 60000 fm/c");
  const Double_t T = 60000.;
  Double_t dP, dL;
  char name[128];

  const Int_t nint = 4;
  const event::integrator ints[nint] = {event::kSymplecticEuler, event::kVelocityVerlet,
                                        event::kRK4, event::kRK45};
  const Double_t dts[nint] = {1., 1., 1., 10.};
  for(Int_t k=0; k<nint; k++){
    event ev; BuildFree(ev);
    ev.SetIntegrator(ints[k]);
    Conservation(ev, TMath::Nint(T/dts[k]), dts[k], dP, dL);
    snprintf(name, sizeof(name), "event %-6s dt=%-3.0f : |dP| / sum|p|", ev.GetIntegratorName().Data(), dts[k]);
    Check(name, dP, 1e-12);
    snprintf(name, sizeof(name), "event %-6s dt=%-3.0f : |dL| / |L|", ev.GetIntegratorName().Data(), dts[k]);
    // central forces : exact (to round-off) for the symplectic schemes only
    if(ints[k]==event::kSymplecticEuler || ints[k]==event::kVelocityVerlet) Check(name, dL, 1e-11);
    else Info(name, dL);
  }
  { // the setup really scatters
    event ev; BuildFree(ev); ev.SetIntegrator("RK4");
    TVector3 v0 = ev.getParticle(0)->GetV();
    ev.DeriveDTN(TMath::Nint(T), 1.);
    Info("deflection of the light particle (deg)", AngleDeg(v0, ev.getParticle(0)->GetV()));
  }
  { // relativistic momentum, gamma m0 v
    event ev; BuildFree(ev); ev.SetIntegrator("RK4"); ev.SetRelativistic(kTRUE);
    Conservation(ev, TMath::Nint(T), 1., dP, dL);
    Check("event RK4 relativistic : |dP| / sum|p|", dP, 1e-12);
  }
  { // pair symmetry off : each force evaluated separately
    event ev; BuildFree(ev); ev.SetIntegrator("RK4"); ev.SetPairSymmetric(kFALSE);
    Conservation(ev, TMath::Nint(T), 1., dP, dL);
    Info("event RK4 pair-symmetric off : |dP| / sum|p|", dP);
  }
  { // original scheme
    event ev; BuildFree(ev);
    Conservation(ev, TMath::Nint(T), 1., dP, dL);
    Info("event TAYLOR dt=1 : |dP| / sum|p|", dP);
  }
  { // sequential update
    oddevent ev; BuildFree(ev);
    Conservation(ev, TMath::Nint(T), 1., dP, dL);
    Info("oddevent dt=1 : |dP| / sum|p|", dP);
    Info("oddevent dt=1 : |dL| / |L|", dL);
  }
}

//------------------------------------------------------------------------------
// 4. Pair symmetry
//------------------------------------------------------------------------------
void BuildMany(event & ev, force * f){
  ev.AddForce(f);
  ev.AddParticle(new EMparticle(206, 79, TVector3(0., 0., 0.), TVector3(0., 0., 0.), true, false));
  ev.AddParticle(new EMparticle(4,  2.3, TVector3(-800., 60., 0.),  TVector3(0.04, 0., 0.), false, false));
  ev.AddParticle(new EMparticle(12, 6.7, TVector3(700., -80., 30.), TVector3(-0.03, 0.002, 0.), false, false));
  ev.AddParticle(new EMparticle(1,  0.9, TVector3(50., 900., -40.), TVector3(0.001, -0.05, 0.), false, false));
  ev.AddParticle(new EMparticle(56, 26, TVector3(3000., 3000., 0.), TVector3(0., 0., 0.), true, false));
}

void CheckPairSymmetry(){
  Header("Pair-symmetric on / off (2 fixed + 3 moving particles, T = 20000 fm/c)");
  const Int_t nint = 4;
  const event::integrator ints[nint] = {event::kSymplecticEuler, event::kVelocityVerlet,
                                        event::kRK4, event::kRK45};
  const Double_t dt = 2.;
  const Long_t   n  = 10000;
  char name[128];
  for(Int_t k=0; k<nint; k++){
    countingCoulomb fon, foff;
    event on;  BuildMany(on, &fon);
    event off; BuildMany(off, &foff);
    on.SetIntegrator(ints[k]);  off.SetIntegrator(ints[k]);
    off.SetPairSymmetric(kFALSE);
    TVector3 xfix0 = on.getParticle(0)->GetX();
    on.DeriveDTN(n, dt);
    off.DeriveDTN(n, dt);
    Double_t dx = 0., dv = 0.;
    Bool_t fixedOK = kTRUE;
    for(Int_t i=0; i<on.getNParticle(); i++){
      dx = TMath::Max(dx, (on.getParticle(i)->GetX() - off.getParticle(i)->GetX()).Mag());
      dv = TMath::Max(dv, (on.getParticle(i)->GetV() - off.getParticle(i)->GetV()).Mag());
    }
    if(on.getParticle(0)->GetX()!=xfix0 || on.getParticle(0)->GetV().Mag()!=0. ||
       off.getParticle(0)->GetX()!=xfix0 || on.getParticle(4)->GetV().Mag()!=0.) fixedOK = kFALSE;
    TString nm = on.GetIntegratorName();
    snprintf(name, sizeof(name), "%-6s : max |x_on - x_off| (fm)", nm.Data());
    Check(name, dx, 1e-8);
    snprintf(name, sizeof(name), "%-6s : max |v_on - v_off| (c)", nm.Data());
    Check(name, dv, 1e-12);
    // 3 moving among 5 : 9 pairs with a moving particle once, 3*4 = 12 ordered pairs
    Double_t ratio = (Double_t) fon.ncall / (Double_t) foff.ncall;
    snprintf(name, sizeof(name), "%-6s : |evaluations on/off - 9/12|  (ratio %.4f)", nm.Data(), ratio);
    Check(name, TMath::Abs(ratio - 9./12.), 1e-12);
    snprintf(name, sizeof(name), "%-6s : invincible particles stay fixed", nm.Data());
    CheckTrue(name, fixedOK);
  }
}

//------------------------------------------------------------------------------
// 5. Ownership
//------------------------------------------------------------------------------
void CheckOwnership(){
  Header("Create / copy / pre-derive / delete events with inspectors");
  const Int_t nloop = 2000;

  event * evT = new event();
  BuildFixed(*evT, 50., -1000., TVector3(0.05, 0., 0.));
  evT->SetIntegrator("VERLET");
  evT->AddInspector(new inspectorP(evT->getParticle(0), evT->getParticle(1), "CNT", 10));
  evT->AddInspector(new inspectorE(evT, "TMG", 20));
  inspectorE * ab = new inspectorE(evT, "DSL", 1.);
  ab->SetAbort();
  evT->AddInspector(ab);

  oddevent * oT = new oddevent();
  BuildFixed(*oT, 50., -1000., TVector3(0.05, 0., 0.));
  oT->AddInspector(new inspectorP(oT->getParticle(0), oT->getParticle(1), "CNT", 10));

  Int_t nok = 0;
  Bool_t linkOK = kTRUE;
  for(Int_t i=0; i<nloop; i++){
    event * ev = new event(*evT);
    if(((inspectorE*) ev->getInspector(1))->GetEvent()!=ev) linkOK = kFALSE;
    if(((inspectorP*) ev->getInspector(0))->GetParticle2()!=ev->getParticle(1)) linkOK = kFALSE;
    if(i%2==0){
      ev->preDerive();
      ev->preDerive();          // replaces the first pre-derived copy
    }
    if(i%4==0) ev->resetPreDerive();
    if(i%3==0){
      event * ev2 = new event(*ev); // copy of a copy
      ev2->DeriveInspect(1);
      delete ev2;
    }
    if(ev->DeriveInspect(1)) nok++;
    ev->AddInspector(new inspectorE(ev, "TMG", 1e9)); // added after the copy : also owned
    delete ev;

    oddevent * oe = new oddevent(*oT);
    oe->DeriveInspect(1);
    delete oe;
  }
  CheckTrue("events copied and deleted without crash", kTRUE,
            TString::Format("%d loops", nloop).Data());
  CheckTrue("copied inspectors watch the copy", linkOK);
  CheckTrue("every copied event ends by its inspectors", nok==nloop,
            TString::Format("%d / %d", nok, nloop).Data());
  CheckTrue("template untouched (no steps, 3 inspectors)",
            evT->GetDeriveN()==0 && evT->getNInspector()==3);
  delete evT;
  delete oT;
  { // stack event with inspectors : deleted at the end of the scope
    event ev; BuildFixed(ev, 50., -1000., TVector3(0.05, 0., 0.));
    ev.AddInspector(new inspectorE(&ev, "TMG", 5));
    ev.DeriveInspect(1);
  }
  CheckTrue("stack event with inspectors goes out of scope", kTRUE);
}

int main(){
  CheckRK45();
  CheckRelativistic();
  CheckConservation();
  CheckPairSymmetry();
  CheckOwnership();

  std::cout<<std::endl<<(gFail==0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")
           <<" ("<<gFail<<" failed)"<<std::endl;
  return gFail==0 ? 0 : 1;
}
