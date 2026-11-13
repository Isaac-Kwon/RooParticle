////////////////////////////////////////////
//
// RooParticle tutorial, part 2 : stopping automatically
//
// The same alpha on a fixed gold-like target as in part 1, but the event is
// derived until its inspectors (terminating conditions) are satisfied.
// Shows the default AND combination, why it is the default, OR mode, abort
// inspectors, the step cap, an undefined method, and GetStopReason().
//
// No arguments. Exit code is non-zero when an event does not end the way
// the tutorial says.
//
////////////////////////////////////////////

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"

#include "TMath.h"
#include "TVector3.h"
#include <cstdio>

static coulombForce gCoulomb;
static Int_t gFail = 0;

// alpha (particle 1) on a fixed target (particle 0) ; impact parameter b [fm]
void Build(event & ev, Double_t b){
  ev.AddForce(&gCoulomb);
  ev.AddParticle(new EMparticle(197., 79., TVector3(0., 0., 0.), TVector3(), true, false));
  ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., b, 0.), TVector3(0.05, 0., 0.), false, false));
  ev.SetIntegrator("VERLET");
}

//! [report]
// One line per event: how it ended, how far the alpha got, and how close its
// speed is to the speed at infinity v_inf = sqrt(2 E / m) (fixed target: all
// the energy E = KE + PE goes to the alpha far away).
void Report(const char* title, event & ev, Bool_t ok){
  particle* target = ev.getParticle(0);
  particle* alpha  = ev.getParticle(1);
  Double_t r    = (alpha->GetX() - target->GetX()).Mag();
  Double_t vinf = TMath::Sqrt(2. * ev.GetNetE() / alpha->GetM0());
  std::printf("%-34s %-9s %8ld steps  r = %7.1f fm  |v|/v_inf = %.5f\n",
              title, ev.GetStopReasonName().Data(), (long) ev.GetDeriveN(), r,
              alpha->GetV().Mag() / vinf);
  //! [report]
  if(!ok){ std::printf("  FAIL : %s\n", title); gFail++; }
}

int main(){
  {
    //! [and]
    event ev;
    Build(ev, 50.);
    particle* target = ev.getParticle(0);
    particle* alpha  = ev.getParticle(1);
    // pair inspectors: angle between r = x_alpha - x_target and v < 1 degree
    // (moving radially outward), and more than 50000 steps of the pair
    ev.AddInspector(new inspectorP(target, alpha, "DEG", 1.));
    ev.AddInspector(new inspectorP(target, alpha, "CNT", 50000));
    Bool_t normal = ev.DeriveInspect();   // = DeriveInspect(1, 1.): inspect every step, dt = 1 fm/c
    Report("b=50: DEG<1 AND CNT>50000", ev, normal && ev.GetStopReason() == event::kStopInspector);
    //! [and]
  }
  {
    //! [why]
    // Nearly head-on: right after the turning point the alpha moves straight
    // out, so "DEG<1" is satisfied at once, close to the target.
    event a;
    Build(a, 1.);
    a.AddInspector(new inspectorP(a.getParticle(0), a.getParticle(1), "DEG", 1.));
    a.DeriveInspect();
    Report("b=1:  DEG<1 alone", a, a.GetDeriveN() < 25000);

    // With the minimum number of steps (AND), it is carried far away.
    event b;
    Build(b, 1.);
    b.AddInspector(new inspectorP(b.getParticle(0), b.getParticle(1), "DEG", 1.));
    b.AddInspector(new inspectorP(b.getParticle(0), b.getParticle(1), "CNT", 50000));
    b.DeriveInspect();
    Report("b=1:  DEG<1 AND CNT>50000", b, b.GetDeriveN() == 50001);
    //! [why]
  }
  {
    //! [event]
    // An event inspector looks at the whole event: here the interaction
    // energy has become small, |PE / KE| < 0.01, AND the pair is receding.
    event ev;
    Build(ev, 50.);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "VRG", 0.));
    ev.AddInspector(new inspectorE(&ev, "PKL", 0.01));
    ev.DeriveInspect();
    Report("b=50: VRG>0 AND PKL<0.01", ev, ev.GetStopReason() == event::kStopInspector);
    //! [event]
  }
  {
    //! [or]
    // OR: stop at whichever comes first.
    event ev;
    Build(ev, 50.);
    ev.SetInspectMode(event::kInspectOR);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", 3000.));
    ev.AddInspector(new inspectorE(&ev, "TMG", 30000.));   // evolved time > 30000 fm/c
    ev.DeriveInspect();
    Report("b=50: DSTG>3000 OR TMG>30000", ev, ev.GetDeriveN() == 30001);
    //! [or]
  }
  {
    //! [abort]
    // An abort inspector ends the event by itself, whatever the others say.
    // DeriveInspect() then returns kFALSE and the stop reason is ABORT.
    event ev;
    Build(ev, 1.);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", 3000.));
    inspector* close = new inspectorE(&ev, "DSL", 60.);   // two particles closer than 60 fm
    close->SetAbort();
    ev.AddInspector(close);
    Bool_t normal = ev.DeriveInspect();
    Report("b=1:  DSTG>3000, abort DSL<60", ev, !normal && ev.GetStopReason() == event::kStopAbort);
    //! [abort]
  }
  {
    //! [cap]
    // The safety cap: at most this many steps per DeriveInspect() call
    // (default 100,000,000; 0 = unlimited). Returns kFALSE, reason MAXSTEPS.
    event ev;
    Build(ev, 50.);
    ev.SetMaxInspectSteps(10000);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DSTG", 3000.));
    Bool_t normal = ev.DeriveInspect();
    Report("b=50: DSTG>3000, cap 10000", ev, !normal && ev.GetStopReason() == event::kStopMaxSteps);
    //! [cap]
  }
  {
    //! [invalid]
    // A misspelt method is refused before anything is derived.
    event ev;
    Build(ev, 50.);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DGE", 1.));
    Bool_t normal = ev.DeriveInspect();
    Report("b=50: \"DGE\" (typo)", ev, !normal && ev.GetStopReason() == event::kStopInvalid
           && ev.GetDeriveN() == 0);
    //! [invalid]
  }
  {
    //! [switch]
    // What to do with each stop reason in a production loop.
    event ev;
    Build(ev, 50.);
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "DEG", 1.));
    ev.AddInspector(new inspectorP(ev.getParticle(0), ev.getParticle(1), "CNT", 50000));
    ev.DeriveInspect(10, 1.);   // inspect every 10 steps (cheaper), dt = 1 fm/c
    switch(ev.GetStopReason()){
      case event::kStopInspector: /* normal end: use the result      */ break;
      case event::kStopAbort:     /* flag or skip this event         */ break;
      case event::kStopMaxSteps:  /* not finished: raise the cap or fix the conditions */ break;
      case event::kStopInvalid:   /* programming error: misspelt method */ break;
      case event::kStopFailed:    /* RK45 gave up (singular force): the state is unusable */ break;
      default: break;
    }
    Report("b=50: same, inspected every 10", ev, ev.GetDeriveN() % 10 == 0
           && ev.GetStopReason() == event::kStopInspector);
    //! [switch]
  }

  std::printf("%s\n", gFail == 0 ? "TUTORIAL 2 OK" : "TUTORIAL 2 FAILED");
  return gFail == 0 ? 0 : 1;
}
