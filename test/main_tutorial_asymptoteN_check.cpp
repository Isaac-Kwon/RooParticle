////////////////////////////////////////////
//
// RooParticle tutorial, part 4b : asymptotic velocities of three and more bodies
//
// The three free bodies of part 4 (12C at rest, an alpha and a proton coming in).
//   A. Stop with the event inspector "ASYD" (remaining deflection below a
//      value, all pairs receding), read the asymptotic velocities with
//      asymptoteN (pairwise), check IsReliable() and the error estimate
//      against the brute-force reference, for three values of ASYD.
//   B. The same accuracy by integrating far and extrapolating (part 4):
//      stop distance and force evaluations needed.
//   C. Incoming asymptotes of the run started at 1000 fm.
//   D. Flags: pairs still approaching, a bound pair, a body bound to a fixed centre.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"
#include "asymptote.hpp"
#include "asymptoteN.hpp"

#include "TMath.h"
#include "TVector3.h"
#include <cstdio>
#include <vector>

static Int_t gFail = 0;
void Expect(const char* what, Bool_t ok){
  if(!ok){ std::printf("  FAIL : %s\n", what); gFail++; }
}

static coulombForce gCoulomb;

// 12C at rest, alpha and proton coming in ; all free (as in part 4 B)
void Build(event & ev){
  ev.AddForce(&gCoulomb);
  ev.AddParticle(new EMparticle(12., 6., TVector3(0., 0., 0.), TVector3(), false, false));
  ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 20., 0.), TVector3(0.05, 0., 0.), false, false));
  ev.AddParticle(new EMparticle(1., 1., TVector3(-900., -60., 10.), TVector3(0.04, 0., 0.), false, false));
  ev.SetIntegrator("RK45");
  ev.SetTolerance(1e-10, 1e-12);
}

// full force evaluations (sweeps over all pairs) an RK45 event has made :
// 6 per substep (accepted or rejected) + 1 per DeriveDT call
Long_t Sweeps(event & ev, Long_t calls){
  return 6 * (ev.GetRK45NSubsteps() + ev.GetRK45NRejected()) + calls;
}

// largest |v_i - ref_i| / |ref_i| over the particles
Double_t MaxRelDiff(const std::vector<TVector3> & v, const std::vector<TVector3> & ref){
  Double_t d = 0.;
  for(size_t i = 0; i < v.size(); i++) d = TMath::Max(d, (v[i] - ref[i]).Mag() / ref[i].Mag());
  return d;
}

std::vector<TVector3> VInf(const asymptoteN & a){
  std::vector<TVector3> v;
  for(Int_t i = 0; i < a.GetNParticle(); i++) v.push_back(a.GetVInf(i));
  return v;
}

std::vector<TVector3> VNow(event & ev){
  std::vector<TVector3> v;
  for(Int_t i = 0; i < ev.getNParticle(); i++) v.push_back(ev.getParticle(i)->GetV());
  return v;
}

int main(){
  // ------------------------------------------------------------------ A
  std::printf("== A. stop with ASYD, then asymptoteN (pairwise) vs brute force ==\n");
  std::printf("%7s %9s %8s %8s %12s %12s %12s %9s %8s\n", "ASYD", "stop (fm)", "steps", "sweeps",
              "true turn", "raw rel err", "pair rel err", "estimate", "reliable");
  const Double_t thr[3] = {1e-2, 1e-3, 1e-4};
  Double_t stopR[3], relPair[3], relRaw[3];
  Long_t   sweeps[3];
  for(Int_t k = 0; k < 3; k++){
    //! [asyd]
    event ev;
    Build(ev);
    ev.AddInspector(new inspectorE(&ev, "ASYD", thr[k])); // every particle turns < thr more (rad), all pairs receding
    inspector* drift = new inspectorE(&ev, "EDG", 1e-6);  // numerical failure
    drift->SetAbort();
    ev.AddInspector(drift);
    ev.SetMaxInspectSteps(1000000);
    ev.DeriveInspect(1, 100.);                           // RK45 subdivides; inspect every 100 fm/c

    asymptoteN pair(&ev);                                // kPairwise, outgoing, from the present state
    asymptoteN brute(&ev, asymptoteN::kBrute);           // reference (private copy, integrated on)
    //! [asyd]
    std::vector<TVector3> vref = VInf(brute), vnow = VNow(ev);
    Double_t turn = 0.;
    for(size_t i = 0; i < vref.size(); i++) turn = TMath::Max(turn, asymptote::Angle(vnow[i], vref[i]));
    stopR[k]   = ev.GetMinDistance();
    sweeps[k]  = Sweeps(ev, ev.GetDeriveN());
    relRaw[k]  = MaxRelDiff(vnow, vref);
    relPair[k] = MaxRelDiff(VInf(pair), vref);
    std::printf("%7.0e %9.0f %8ld %8ld %12.2e %12.1e %12.1e %9.1e %8s\n", thr[k], stopR[k],
                (long) ev.GetDeriveN(), sweeps[k], turn, relRaw[k], relPair[k], pair.GetMaxRelError(),
                pair.IsReliable() ? "yes" : "no");
    Expect("A ends by ASYD", ev.GetStopReason() == event::kStopInspector);
    Expect("A reliable", pair.IsReliable());
    Expect("A estimate covers the actual error", relPair[k] < pair.GetMaxRelError());
    Expect("A true remaining deflection is of the ASYD size", turn < 1.5 * thr[k] && turn > 0.5 * thr[k]);

    if(k == 0){
      std::printf("\ndetails at ASYD = %.0e:\n", thr[k]);
      //! [pairwise]
      std::printf("stop                : %s, closest pair %.0f fm, remaining deflection (first order) %.3e rad\n",
                  ev.GetStopReasonName().Data(), ev.GetMinDistance(), asymptoteN::RemainingDeflection(&ev));
      std::printf("pairwise            : reliable %s, flags [%s], largest relative error %.1e\n",
                  pair.IsReliable() ? "yes" : "no", pair.GetFlagNames(pair.GetFlags()).Data(), pair.GetMaxRelError());
      for(Int_t i = 0; i < ev.getNParticle(); i++){
        std::printf("  particle %d : |v| now %.9f c, |v_inf| %.9f c +- %.1e, still to turn %.2e rad,"
                    " pairwise - brute %.1e c\n",
                    i, ev.getParticle(i)->GetV().Mag(), pair.GetVInf(i).Mag(), pair.GetError(i),
                    pair.GetDeflection(i), (pair.GetVInf(i) - brute.GetVInf(i)).Mag());
      }
      std::printf("brute force         : %ld RK45 substeps (%ld sweeps), integrated %.1e fm/c on\n",
                  (long) brute.GetNSubsteps(), (long) brute.GetNForceSweeps(), brute.GetIntegratedTime());
      //! [pairwise]

      // ---------------------------------------------------------------- B
      std::printf("\n== B. the same accuracy by integrating far from this stop ==\n");
      //! [far]
      // Continue a copy with the event's own RK45 to larger and larger distance;
      // compare the raw velocities and the 1/R (Richardson) extrapolation with
      // the brute-force reference. Sweeps are counted from the ASYD stop on.
      event far(ev);
      std::printf("pairwise at the stop: relative error %.1e, 0 force sweeps (closed form)\n", relPair[k]);
      std::printf("%9s %10s %14s %14s\n", "R (fm)", "sweeps", "raw rel err", "Richardson");
      std::vector<TVector3> prev;
      Long_t calls = 0;
      Double_t Rrich = 0., Rraw = 0.;
      Long_t   Srich = 0,  Sraw = 0;
      for(Double_t R = 1e4; R <= 1e9 * 1.0001; R *= 10.){
        while(far.GetMinDistance() < R){ far.DeriveDT(R / 100.); calls++; }
        std::vector<TVector3> v = VNow(far);
        Double_t eraw = MaxRelDiff(v, vref), erich = -1.;
        if(!prev.empty()){
          std::vector<TVector3> vr;
          for(size_t i = 0; i < v.size(); i++) vr.push_back(v[i] + (v[i] - prev[i]) * (1. / 9.));
          erich = MaxRelDiff(vr, vref);
        }
        std::printf("%9.0e %10ld %14.1e %14.1e\n", R, Sweeps(far, calls), eraw, erich);
        if(Rrich == 0. && erich >= 0. && erich <= relPair[k]){ Rrich = R; Srich = Sweeps(far, calls); }
        if(Rraw  == 0. && eraw <= relPair[k])              { Rraw  = R; Sraw  = Sweeps(far, calls); }
        prev = v;
      }
      //! [far]
      std::printf("pairwise accuracy reached: extrapolated at R = %.0e fm after %ld more sweeps,"
                  " raw at R = %.0e fm after %ld more sweeps\n", Rrich, Srich, Rraw, Sraw);
      Expect("B extrapolation needs R >= 1e5 fm", Rrich >= 1e5);
      Expect("B raw needs R >= 1e6 fm", Rraw >= 1e6);
    }
  }
  Expect("A smaller ASYD : farther stop, smaller raw error", stopR[2] > stopR[1] && stopR[1] > stopR[0] && relRaw[2] < relRaw[0]);
  Expect("A pairwise error falls like ASYD^2", relPair[0] / relPair[2] > 3e3 && relPair[0] / relPair[2] < 3e4);

  {
    //! [asyv]
    // ASYV : the relative change of speed and direction, |v_inf - v| / |v_inf|, still to come
    event ev;
    Build(ev);
    ev.AddInspector(new inspectorE(&ev, "ASYV", 1e-3));
    ev.DeriveInspect(1, 100.);
    asymptoteN brute(&ev, asymptoteN::kBrute);
    Double_t change = MaxRelDiff(VNow(ev), VInf(brute));
    std::printf("\nASYV < 1e-3         : %s at %.0f fm, true remaining relative change %.2e\n",
                ev.GetStopReasonName().Data(), ev.GetMinDistance(), change);
    //! [asyv]
    Expect("ASYV true change of its size", change < 1.5e-3 && change > 0.5e-3);
  }

  // ------------------------------------------------------------------ C
  std::printf("\n== C. incoming asymptotes of the run started at 1000 fm ==\n");
  //! [incoming]
  event start;
  Build(start);                                       // the initial state only, not derived
  asymptoteN in(&start, asymptoteN::kPairwise, asymptoteN::kIncoming);
  asymptoteN inRef(&start, asymptoteN::kBrute, asymptoteN::kIncoming);
  for(Int_t i = 1; i < start.getNParticle(); i++){
    TVector3 v0 = start.getParticle(i)->GetV();
    std::printf("  particle %d : nominal |v| %.6f c, incoming |v| %.9f c, tilt %.4f deg,"
                " pairwise - brute %.1e c (estimate %.1e)\n",
                i, v0.Mag(), in.GetVInf(i).Mag(), asymptote::Angle(v0, in.GetVInf(i)) * TMath::RadToDeg(),
                (in.GetVInf(i) - inRef.GetVInf(i)).Mag(), in.GetError(i));
  }
  std::printf("  12C        : at rest at the start, incoming |v| %.2e c\n", in.GetVInf(0).Mag());
  //! [incoming]
  for(Int_t i = 0; i < 3; i++)
    Expect("C incoming within the estimate", (in.GetVInf(i) - inRef.GetVInf(i)).Mag() < in.GetError(i) + inRef.GetError(i));

  // ------------------------------------------------------------------ D
  std::printf("\n== D. flags ==\n");
  //! [flags]
  { // the initial state, outgoing: every pair is still approaching
    asymptoteN a(&start);
    std::printf("outgoing from t = 0 : separating %s, flags [%s], reliable %s, largest error %.1e c\n",
                a.IsSeparating() ? "yes" : "no", a.GetFlagNames(a.GetFlags()).Data(),
                a.IsReliable() ? "yes" : "no", a.GetMaxError());
    Expect("D approaching flagged", (a.GetFlags() & asymptoteN::kApproaching) && !a.IsSeparating());
  }
  { // a bound pair (+1, -1, 10 fm apart) drifting at 0.02 c, an alpha passing at 200 fm
    event b;
    b.AddForce(&gCoulomb);
    b.AddParticle(new EMparticle(1.,  1., TVector3(0.,  5., 0.), TVector3(0.026, 0., 0.), false, false));
    b.AddParticle(new EMparticle(1., -1., TVector3(0., -5., 0.), TVector3(0.014, 0., 0.), false, false));
    b.AddParticle(new EMparticle(4.,  2., TVector3(1000., 200., 0.), TVector3(-0.05, 0., 0.), false, false));
    asymptoteN a(&b);
    std::printf("bound pair          : flags [%s] [%s] [%s], clusters %d, cluster |V_inf| %.6f c, reliable %s\n",
                a.GetFlagNames(a.GetFlags(0)).Data(), a.GetFlagNames(a.GetFlags(1)).Data(),
                a.GetFlagNames(a.GetFlags(2)).Data(), a.GetNClusters(),
                a.GetNClusters() > 0 ? a.GetClusterVInf(0).Mag() : 0., a.IsReliable() ? "yes" : "no");
    Expect("D cluster", a.GetNClusters() == 1 && (a.GetFlags(0) & asymptoteN::kInCluster));
  }
  { // a light negative particle bound to a fixed Au, an alpha passing
    event f;
    f.AddForce(&gCoulomb);
    f.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
    f.AddParticle(new EMparticle(1., -1., TVector3(100., 0., 0.), TVector3(0., 0.02, 0.), false, false));
    f.AddParticle(new EMparticle(4., 2., TVector3(-1000., 300., 0.), TVector3(0.05, 0., 0.), false, false));
    asymptoteN a(&f);
    std::printf("bound to fixed      : light particle [%s], v_inf %.1f c ; alpha [%s]\n",
                a.GetFlagNames(a.GetFlags(1)).Data(), a.GetVInf(1).Mag(), a.GetFlagNames(a.GetFlags(2)).Data());
    Expect("D bound to fixed", (a.GetFlags(1) & asymptoteN::kBoundToFixed) && a.GetVInf(1).Mag() == 0.);
  }
  //! [flags]

  std::printf("%s\n", gFail == 0 ? "TUTORIAL 4b OK" : "TUTORIAL 4b FAILED");
  return gFail == 0 ? 0 : 1;
}
