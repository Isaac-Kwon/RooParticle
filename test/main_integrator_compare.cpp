////////////////////////////////////////////
//
// RooParticle Test Application
// Integrator and terminating condition check
// (Fixed Target, one alpha particle on Au)
//
// 1. Every integrator is run with several dt on the same event.
//    Energy drift and error of the final velocity direction are printed,
//    with the observed order of convergence.
// 2. Every kind of terminating condition is run once.
//
// Exit code is non-zero when a check fails.
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

// Target (invincible) at origin and one incident alpha particle.
// The event owns the particles.
void BuildEvent(event & ev, Double_t b, Double_t x0 = -1000.){
  ev.AddForce(&gCoulomb);
  ev.AddParticle(new EMparticle(206, 79, TVector3(0.,0.,0.), TVector3(0.,0.,0.), true, false));
  ev.AddParticle(new EMparticle(4, 2, TVector3(x0, b, 0.), TVector3(0.05, 0., 0.), false, false));
}

struct runResult{
  Double_t dE;   // relative energy drift
  TVector3 v;    // final velocity of the incident particle
};

runResult Run(event::integrator integ, Double_t dt, Double_t T, Double_t b){
  event ev;
  BuildEvent(ev, b);
  ev.SetIntegrator(integ);

  Double_t e0 = ev.GetNetE();
  ev.DeriveDTN(TMath::Nint(T/dt), dt);
  Double_t e1 = ev.GetNetE();

  runResult r;
  r.dE = TMath::Abs((e1-e0)/e0);
  r.v  = ev.getParticle(1)->GetV();
  return r;
}

Double_t AngleDeg(TVector3 a, TVector3 b){
  return TMath::ATan2(a.Cross(b).Mag(), a.Dot(b)) * TMath::RadToDeg();
}

Int_t CheckIntegrators(){
  const Double_t b = 50.;      // impact parameter (fm)
  const Double_t T = 40000.;   // total time (fm/c)
  const Int_t    ndt = 3;
  // coarse enough that the error of RK4 stays above the round-off level
  const Double_t dts[ndt] = {32., 16., 8.};

  const Int_t nint = 4;
  const event::integrator ints[nint] = {event::kTaylor, event::kSymplecticEuler,
                                        event::kVelocityVerlet, event::kRK4};
  const Double_t expected[nint] = {1., 1., 2., 4.};

  // reference : RK4 with a small step
  runResult ref = Run(event::kRK4, 0.125, T, b);
  TVector3 vin(0.05, 0., 0.);

  std::cout<<"== Integrators (b = "<<b<<" fm, T = "<<T<<" fm/c) =="<<std::endl;
  std::cout<<"reference scattering angle (RK4, dt=0.125) : "<<AngleDeg(vin, ref.v)<<" deg"<<std::endl;
  printf("%-8s %6s %14s %16s %8s\n", "method", "dt", "|dE/E|", "angle err (deg)", "order");

  Int_t nfail = 0;
  for(Int_t i=0; i<nint; i++){
    event ev;
    ev.SetIntegrator(ints[i]);
    TString name = ev.GetIntegratorName();

    Double_t preverr = -1.;
    Double_t order = 0.;
    for(Int_t j=0; j<ndt; j++){
      runResult r = Run(ints[i], dts[j], T, b);
      Double_t err = AngleDeg(ref.v, r.v);
      if(preverr>0. && err>0.){
        order = TMath::Log2(preverr/err);
        printf("%-8s %6.3f %14.3e %16.3e %8.2f\n", name.Data(), dts[j], r.dE, err, order);
      }else{
        printf("%-8s %6.3f %14.3e %16.3e %8s\n", name.Data(), dts[j], r.dE, err, "-");
      }
      preverr = err;
    }
    // observed order (from the two smallest dt) should match the scheme
    if(TMath::Abs(order-expected[i])>0.5){
      std::cout<<"  FAIL : expected order "<<expected[i]<<std::endl;
      nfail++;
    }
  }
  return nfail;
}

// One terminating setup : report how DeriveInspect() ended.
Int_t Report(const char* title, event & ev, Bool_t rtn, event::stopReason expected){
  particle* p = ev.getParticle(1);
  Double_t r = (p->GetX() - ev.getParticle(0)->GetX()).Mag();
  Bool_t ok = (ev.GetStopReason()==expected);
  printf("%-44s %-10s steps %8ld  time %9.1f  r %9.1f  %s\n",
         title, ev.GetStopReasonName().Data(), (long) ev.GetDeriveN(), ev.GetTime(), r,
         ok ? "ok" : "FAIL");
  if(rtn != (expected==event::kStopInspector)) ok = kFALSE;
  return ok ? 0 : 1;
}

Int_t CheckTermination(){
  const Double_t b = 50.;
  Int_t nfail = 0;

  std::cout<<std::endl<<"== Terminating conditions (VERLET, dt = 1) =="<<std::endl;

  { // original setup : asymptote angle AND minimum number of steps
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DEG", 1));
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "CNT", 50000));
    nfail += Report("DEG<1 AND CNT>50000", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // same inspectors, OR mode
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.SetInspectMode(event::kInspectOR);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DEG", 1));
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "CNT", 50000));
    nfail += Report("DEG<1 OR CNT>50000", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // receding and far away
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "VRG", 0));
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", 3000));
    nfail += Report("receding (VRG>0) AND DSTG>3000", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // receding and interaction energy negligible
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "VRG", 0));
    ev.AddInspector(new inspectorE(&ev, "PKL", 0.01));
    nfail += Report("receding AND |PE/KE|<0.01 (PKL)", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // receding and force negligible
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "VRG", 0));
    ev.AddInspector(new inspectorE(&ev, "ACL", 1e-9));
    nfail += Report("receding AND |a|<1e-9 (ACL)", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // evolved time
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.AddInspector(new inspectorE(&ev, "TMG", 30000));
    nfail += Report("time>30000 (TMG)", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // leaving the event volume
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    eventVGeneral vol(4000., 4000., 4000.);
    ev.SetVolume(&vol);
    ev.AddInspector(new inspectorE(&ev, "VOUT"));
    nfail += Report("out of 4000 fm box (VOUT)", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // abort inspector : close encounter stops the event before the others
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", 3000));
    inspectorE* ab = new inspectorE(&ev, "DSL", 200);
    ab->SetAbort();
    ev.AddInspector(ab);
    nfail += Report("DSTG>3000, abort on DSL<200", ev, ev.DeriveInspect(1), event::kStopAbort);
  }
  { // abort inspector : energy drift (coarse step with the original scheme)
    event ev; BuildEvent(ev, 1.); ev.SetIntegrator("TAYLOR");
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", 3000));
    inspectorE* ab = new inspectorE(&ev, "EDG", 1e-4);
    ab->SetAbort();
    ev.AddInspector(ab);
    nfail += Report("TAYLOR dt=20, abort on drift>1e-4 (EDG)", ev, ev.DeriveInspect(1, 20.), event::kStopAbort);
  }
  { // safety cap
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.SetMaxInspectSteps(10000);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", 3000));
    nfail += Report("DSTG>3000, max steps 10000", ev, ev.DeriveInspect(1), event::kStopMaxSteps);
  }
  { // undefined method
    event ev; BuildEvent(ev, b); ev.SetIntegrator("VERLET");
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DGE", 1));
    nfail += Report("undefined method \"DGE\"", ev, ev.DeriveInspect(1), event::kStopInvalid);
  }
  { // oddevent with an invincible target and a CNT inspector
    oddevent ev; BuildEvent(ev, b);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DEG", 1));
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "CNT", 50000));
    nfail += Report("oddevent : DEG<1 AND CNT>50000", ev, ev.DeriveInspect(1), event::kStopInspector);
  }
  { // copied event keeps the integrator and inspects itself
    event evT; BuildEvent(evT, b); evT.SetIntegrator("RK4");
    evT.AddInspector(new inspectorE(&evT, "TMG", 30000));
    event ev(evT);
    Bool_t rtn = ev.DeriveInspect(1);
    Int_t f = Report("copy of a template event (RK4, TMG)", ev, rtn, event::kStopInspector);
    if(ev.GetIntegrator()!=event::kRK4 || evT.GetDeriveN()!=0){
      std::cout<<"  FAIL : copied event"<<std::endl;
      f = 1;
    }
    nfail += f;
  }
  return nfail;
}

int main(){
  Int_t nfail = 0;
  nfail += CheckIntegrators();
  nfail += CheckTermination();

  std::cout<<std::endl<<(nfail==0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")<<std::endl;
  return nfail==0 ? 0 : 1;
}
