////////////////////////////////////////////
//
// RooParticle Test Application
// Analytic two-body asymptote (asymptote.hpp) against RK4 trajectories
//
// "Truth" : RK4, dt = 0.5 fm/c while the pair is closer than 5000 fm, then
// dt = 1e-3 r / |v| out to r = 1e14 fm (the residual deflection there is
// ~ a / r ~ 1e-13 rad).
//
// 1. invariance   : outgoing direction, v_inf and effective b extrapolated
//                   from many instants of one trajectory agree.
// 2. consistency  : extrapolated angle = Rutherford formula at (b_eff, v_inf),
//                   closest approach = formula = integrated minimum.
// 3. benefit      : extrapolating from a short stop (300 fm, 1000 fm after
//                   closest approach) matches the long integration; the raw
//                   direction there does not.
// 4. free target  : lab velocities of both bodies, momentum and energy.
// 5. head-on and attractive cases, and the status of special states.
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

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"

// Coulomb constant (amu fm c^2 / e^2). Same number as coulombForce::coulombK,
// which has no getter.
static const Double_t kCoulombK = 0.00154629;

static coulombForce gCoulomb;
static Int_t gFail = 0;
static Int_t gCheck = 0;

static const Double_t kRNear  = 5000.;  // fixed step inside this distance
static const Double_t kDtNear = 0.5;
static const Double_t kFrac   = 1e-3;   // dt = kFrac r/|v| outside
static const Double_t kRFar   = 1e14;   // end of the long integration

// error measure "val" must be below "tol"
void Check(const char* what, Double_t val, Double_t tol){
  Bool_t ok = TMath::Finite(val) && val <= tol;
  printf("  %-62s %12.3e < %8.1e  %s\n", what, val, tol, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
// "val" must be above "low" (e.g. the raw error is really there)
void CheckAbove(const char* what, Double_t val, Double_t low){
  Bool_t ok = TMath::Finite(val) && val >= low;
  printf("  %-62s %12.3e > %8.1e  %s\n", what, val, low, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
void CheckTrue(const char* what, Bool_t ok){
  printf("  %-62s %25s  %s\n", what, "", ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
void Info(const char* what, Double_t val, const char* unit=""){
  printf("  %-62s %18.10g %s\n", what, val, unit);
}

Double_t RelDiff(Double_t a, Double_t b){ return TMath::Abs(a-b)/TMath::Abs(b); }
Double_t RelDiff(const TVector3 & a, const TVector3 & b, Double_t scale){ return (a-b).Mag()/scale; }

// projectile is particle 0, target particle 1
event* Build(Double_t m1, Double_t q1, Double_t m2, Double_t q2, Bool_t fixed2,
             Double_t b, Double_t v0, Double_t x0 = -1000.){
  event* ev = new event();
  ev->AddForce(&gCoulomb);
  ev->AddParticle(new EMparticle(m1, q1, TVector3(x0, b, 0.), TVector3(v0, 0., 0.), false, false));
  ev->AddParticle(new EMparticle(m2, q2, TVector3(0., 0., 0.), TVector3(0., 0., 0.), fixed2, false));
  ev->SetIntegrator("RK4");
  return ev;
}

struct snap{
  Double_t t;
  TVector3 x1, v1, x2, v2;
};

snap Take(event* ev, Double_t t){
  snap s;
  s.t = t;
  s.x1 = ev->getParticle(0)->GetX(); s.v1 = ev->getParticle(0)->GetV();
  s.x2 = ev->getParticle(1)->GetX(); s.v2 = ev->getParticle(1)->GetV();
  return s;
}

asymptote FromSnap(const snap & s, event* ev, Double_t kappa){
  particle* p1 = ev->getParticle(0);
  particle* p2 = ev->getParticle(1);
  return asymptote(s.x1, s.v1, p1->GetM0(), p1->IsInvincible(),
                   s.x2, s.v2, p2->GetM0(), p2->IsInvincible(), kappa);
}

// Integrate until the pair is receding and farther than rstop.
// rmin : smallest distance so far, refined by a parabola through three steps.
// snaps : if given, the state every "every" steps of the fixed-step phase.
struct runState{
  Double_t t    = 0.;
  Double_t rmin = 1e300;
  Double_t r[3] = {0., 0., 0.};  // last three distances of the fixed-step phase
  Int_t    nr   = 0;
  Long_t   n    = 0;
};

void Integrate(event* ev, runState & rs, Double_t rstop,
               std::vector<snap>* snaps = nullptr, Long_t every = 1000){
  particle* p1 = ev->getParticle(0);
  particle* p2 = ev->getParticle(1);
  for(Long_t guard=0; guard<100000000; guard++){
    TVector3 dx = p1->GetX() - p2->GetX();
    TVector3 dv = p1->GetV() - p2->GetV();
    Double_t r  = dx.Mag();
    if(r < rs.rmin) rs.rmin = r;
    if(dx.Dot(dv) > 0. && r >= rstop) return;

    Double_t dt;
    if(r < kRNear){
      dt = kDtNear;
      if(snaps && rs.n % every == 0) snaps->push_back(Take(ev, rs.t));
      rs.r[0] = rs.r[1]; rs.r[1] = rs.r[2]; rs.r[2] = r; rs.nr++;
      if(rs.nr >= 3 && rs.r[1] <= rs.r[0] && rs.r[1] <= rs.r[2]){
        Double_t den = rs.r[0] - 2.*rs.r[1] + rs.r[2];
        if(den > 0.){
          Double_t rv = rs.r[1] - (rs.r[2]-rs.r[0])*(rs.r[2]-rs.r[0])/(8.*den);
          if(rv < rs.rmin) rs.rmin = rv;
        }
      }
    }else{
      dt = kFrac * r / dv.Mag();
      rs.nr = 0;
    }
    ev->DeriveDTN(1, dt);
    rs.t += dt;
    rs.n++;
  }
  std::cout<<"  integration did not end"<<std::endl;
  gFail++;
}

////////////////////////////////////////////
// 1. invariance along one trajectory
////////////////////////////////////////////
void Invariance(){
  const Double_t m1 = 4., q1 = 2., m2 = 206., q2 = 79.;
  const Double_t b0 = 50., v0 = 0.05;
  const Double_t kappa = kCoulombK*q1*q2;

  std::cout<<std::endl<<"== 1. Invariance along one trajectory (alpha on fixed Au, b = "
           <<b0<<" fm, v = "<<v0<<" c, from x = -1000 fm) =="<<std::endl;

  event* ev = Build(m1, q1, m2, q2, kTRUE, b0, v0);
  runState rs;
  std::vector<snap> snaps;
  Integrate(ev, rs, 4000., &snaps, 1000);

  asymptote ref = FromSnap(snaps[0], ev, kappa);
  // the sample nearest to the closest approach
  size_t imin = 0;
  for(size_t i=0; i<snaps.size(); i++){
    if((snaps[i].x1-snaps[i].x2).Mag() < (snaps[imin].x1-snaps[imin].x2).Mag()) imin = i;
  }

  printf("  %8s %9s %9s %4s %11s %11s %11s %11s %16s %13s\n",
         "t(fm/c)", "r(fm)", "|PE/KE|", "dir", "d(out)rad", "d(in)rad", "dvinf/vinf", "db/b", "thetaCM(deg)", "rmin(fm)");
  Double_t maxOut = 0., maxIn = 0., maxV = 0., maxB = 0., maxE = 0., maxL = 0., maxA = 0., maxRmin = 0.;
  for(size_t i=0; i<snaps.size(); i++){
    asymptote a = FromSnap(snaps[i], ev, kappa);
    Double_t r   = (snaps[i].x1-snaps[i].x2).Mag();
    Double_t ke  = 0.5*m1*snaps[i].v1.Mag2();
    Double_t dO  = asymptote::Angle(a.GetVOutRel(), ref.GetVOutRel());
    Double_t dI  = asymptote::Angle(a.GetVInRel(),  ref.GetVInRel());
    Double_t dV  = RelDiff(a.GetVinf(), ref.GetVinf());
    Double_t dB  = RelDiff(a.GetImpactParameter(), ref.GetImpactParameter());
    maxOut = TMath::Max(maxOut, dO);
    maxIn  = TMath::Max(maxIn,  dI);
    maxV   = TMath::Max(maxV,   dV);
    maxB   = TMath::Max(maxB,   dB);
    maxE   = TMath::Max(maxE,   RelDiff(a.GetEnergy(), ref.GetEnergy()));
    maxL   = TMath::Max(maxL,   RelDiff(a.GetL(), ref.GetL(), ref.GetLMag()));
    maxA   = TMath::Max(maxA,   RelDiff(a.GetLRL(), ref.GetLRL(), ref.GetLRL().Mag()));
    maxRmin= TMath::Max(maxRmin,RelDiff(a.GetClosestApproach(), ref.GetClosestApproach()));
    if(a.GetStatus()!=asymptote::kUnbound){ gFail++; std::cout<<"  FAIL status "<<a.GetStatusName()<<std::endl; }
    if(i%12==0 || i==imin || i+1==snaps.size()){
      printf("  %8.0f %9.2f %9.2e %4s %11.2e %11.2e %11.2e %11.2e %16.11f %13.9f\n",
             snaps[i].t, r, TMath::Abs(kappa/r)/ke, a.IsOutgoing() ? "out" : "in",
             dO, dI, dV, dB, a.GetThetaCM()*TMath::RadToDeg(), a.GetClosestApproach());
    }
  }
  printf("  (%d instants checked, every 500 fm/c, every 12th printed; deviations from the instant t = 0)\n",
         (Int_t) snaps.size());
  Check("max angle between extrapolated outgoing directions (rad)", maxOut, 1e-10);
  Check("max angle between extrapolated incoming directions (rad)", maxIn, 1e-10);
  Check("max relative deviation of v_inf", maxV, 1e-11);
  Check("max relative deviation of effective b", maxB, 1e-11);
  Check("max relative deviation of E", maxE, 1e-11);
  Check("max relative deviation of L vector", maxL, 1e-11);
  Check("max relative deviation of LRL vector", maxA, 1e-10);
  Check("max relative deviation of closest approach", maxRmin, 1e-11);
  Check("closest approach vs integrated minimum (relative)",
        RelDiff(rs.rmin, ref.GetClosestApproach()), 1e-7);

  std::cout<<"  nominal vs asymptotic incoming state (finite start at 1000 fm):"<<std::endl;
  Info("nominal b (fm)", b0);
  Info("effective b = |L|/(mu v_inf) (fm)", ref.GetImpactParameter());
  Info("nominal v (c)", v0);
  Info("v_inf (c)", ref.GetVinf());
  Info("tilt of the incoming asymptote from the x axis (deg)",
       asymptote::Angle(ref.GetVInRel(), TVector3(1.,0.,0.))*TMath::RadToDeg());
  delete ev;
}

////////////////////////////////////////////
// 2, 3. consistency and benefit, fixed target
////////////////////////////////////////////
void FixedTarget(const char* title, Double_t m1, Double_t q1, Double_t m2, Double_t q2,
                 Double_t b0, Double_t v0){
  const Double_t kappa = kCoulombK*q1*q2;
  std::cout<<std::endl<<"== "<<title<<" (m "<<m1<<", q "<<q1<<" on fixed m "<<m2<<", q "<<q2
           <<"; b = "<<b0<<" fm, v = "<<v0<<" c) =="<<std::endl;

  event* ev = Build(m1, q1, m2, q2, kTRUE, b0, v0);
  runState rs;
  asymptote a0(ev->getParticle(0), ev->getParticle(1), kappa);
  Integrate(ev, rs, 300.);
  snap s300 = Take(ev, rs.t);
  asymptote a300 = FromSnap(s300, ev, kappa);
  Integrate(ev, rs, 1000.);
  snap s1000 = Take(ev, rs.t);
  asymptote a1000 = FromSnap(s1000, ev, kappa);
  Integrate(ev, rs, kRFar);
  TVector3 vfar = ev->getParticle(0)->GetV();
  Double_t rfar = (ev->getParticle(0)->GetX() - ev->getParticle(1)->GetX()).Mag();
  delete ev;

  const Double_t mu = a0.GetMu();
  printf("  %-62s %18s\n", "status", a0.GetStatusName().Data());
  Info("effective b (fm)", a0.GetImpactParameter());
  Info("v_inf (c)", a0.GetVinf());
  Info("eccentricity", a0.GetEccentricity());
  Info("extrapolated thetaCM (deg)", a0.GetThetaCM()*TMath::RadToDeg());
  Info("Rutherford at nominal (b, v) (deg)",
       TMath::Abs(asymptote::RutherfordAngle(b0, v0, mu, kappa))*TMath::RadToDeg());
  Info("long integration angle (deg), r = 1e14 fm",
       asymptote::Angle(TVector3(1.,0.,0.), vfar)*TMath::RadToDeg());
  Info("closest approach (fm)", a0.GetClosestApproach());
  Info("long integration final distance (fm)", rfar);

  // 2. consistency
  Double_t thVec = asymptote::Angle(a0.GetVInRel(), a0.GetVOutRel());
  Double_t thRuth = TMath::Abs(asymptote::RutherfordAngle(a0.GetImpactParameter(), a0.GetVinf(), mu, kappa));
  Check("|theta(in,out vectors) - thetaCM| (rad)", TMath::Abs(thVec - a0.GetThetaCM()), 1e-13);
  Check("|thetaCM - Rutherford(b_eff, v_inf)| (rad)", TMath::Abs(a0.GetThetaCM() - thRuth), 1e-13);
  if(b0 > 0.){
    Check("Rutherford b(thetaCM) vs effective b (relative)",
          RelDiff(asymptote::RutherfordImpactParameter(a0.GetThetaCM(), a0.GetVinf(), mu, kappa),
                  a0.GetImpactParameter()), 1e-12);
  }
  Check("closest approach vs formula at (b_eff, v_inf) (relative)",
        RelDiff(a0.GetClosestApproach(),
                asymptote::RutherfordClosestApproach(a0.GetImpactParameter(), a0.GetVinf(), mu, kappa)), 1e-12);
  Check("closest approach vs integrated minimum (relative)",
        RelDiff(rs.rmin, a0.GetClosestApproach()), 1e-7);
  Check("|v_inf - |v| at r = 1e14 fm| / v_inf", RelDiff(vfar.Mag(), a0.GetVinf()), 1e-11);
  Info("nominal-formula bias |Rutherford(b0,v0) - thetaCM| (deg)",
       TMath::Abs(TMath::Abs(asymptote::RutherfordAngle(b0, v0, mu, kappa)) - a0.GetThetaCM())*TMath::RadToDeg());

  // 3. benefit
  Double_t ext0    = asymptote::Angle(a0.GetVOutRel(),    vfar);
  Double_t ext300  = asymptote::Angle(a300.GetVOutRel(),  vfar);
  Double_t ext1000 = asymptote::Angle(a1000.GetVOutRel(), vfar);
  Double_t raw300  = asymptote::Angle(s300.v1,  vfar);
  Double_t raw1000 = asymptote::Angle(s1000.v1, vfar);
  Check("extrapolated from t = 0 vs long direction (rad)", ext0, 1e-10);
  Check("extrapolated from r = 300 fm (out) vs long direction (rad)", ext300, 1e-10);
  Check("extrapolated from r = 1000 fm (out) vs long direction (rad)", ext1000, 1e-10);
  Check("extrapolated v_inf at r = 300 fm vs long |v| (relative)", RelDiff(a300.GetVinf(), vfar.Mag()), 1e-11);
  printf("  raw direction at r = 300 fm  : off by %.3e rad = %.4f deg, |v| off by %.3e (relative)\n",
         raw300, raw300*TMath::RadToDeg(), RelDiff(s300.v1.Mag(), vfar.Mag()));
  printf("  raw direction at r = 1000 fm : off by %.3e rad = %.4f deg, |v| off by %.3e (relative)\n",
         raw1000, raw1000*TMath::RadToDeg(), RelDiff(s1000.v1.Mag(), vfar.Mag()));
  if(b0 > 0.){
    CheckAbove("raw error / extrapolated error at r = 1000 fm", raw1000/TMath::Max(ext1000, 1e-16), 1e3);
  }
}

////////////////////////////////////////////
// 4. free target : lab velocities, momentum, energy
////////////////////////////////////////////
void FreeTarget(const char* title, Double_t m1, Double_t q1, Double_t m2, Double_t q2,
                Double_t b0, Double_t v0, Bool_t angles){
  const Double_t kappa = kCoulombK*q1*q2;
  std::cout<<std::endl<<"== "<<title<<" (m "<<m1<<", q "<<q1<<" on free m "<<m2<<", q "<<q2
           <<"; b = "<<b0<<" fm, v = "<<v0<<" c) =="<<std::endl;

  event* ev = Build(m1, q1, m2, q2, kFALSE, b0, v0);
  particle* p1 = ev->getParticle(0);
  particle* p2 = ev->getParticle(1);
  // Galilean shift so that the target is at rest at t -> -inf (not at t = 0)
  TVector3 shift = asymptote(p1, p2, kappa).GetV2In();
  p1->applyV(p1->GetV() - shift);
  p2->applyV(p2->GetV() - shift);
  Info("target speed at t = 0 after the shift (c)", shift.Mag());

  asymptote a0(p1, p2, kappa);
  TVector3 P0 = p1->GetV()*m1 + p2->GetV()*m2;
  Double_t E0 = ev->GetNetE();
  Info("|v2 in| asymptotic (c)", a0.GetV2In().Mag());

  runState rs;
  Integrate(ev, rs, 300.);
  asymptote a300(p1, p2, kappa);
  TVector3 raw1 = p1->GetV(), raw2 = p2->GetV();
  Integrate(ev, rs, kRFar);
  TVector3 f1 = p1->GetV(), f2 = p2->GetV();
  TVector3 Pf = f1*m1 + f2*m2;
  delete ev;

  Info("thetaCM (deg)", a0.GetThetaCM()*TMath::RadToDeg());
  Info("thetaLab1 (deg)", a0.GetThetaLab1()*TMath::RadToDeg());
  Info("thetaLab2 (deg)", a0.GetThetaLab2()*TMath::RadToDeg());
  Info("|v1 out| (c)", a0.GetV1Out().Mag());
  Info("|v2 out| (c)", a0.GetV2Out().Mag());

  Check("|v1out(t=0) - v1(long)| / v0", RelDiff(a0.GetV1Out(), f1, v0), 1e-10);
  Check("|v2out(t=0) - v2(long)| / v0", RelDiff(a0.GetV2Out(), f2, v0), 1e-10);
  Check("|v1out(r=300) - v1(long)| / v0", RelDiff(a300.GetV1Out(), f1, v0), 1e-10);
  Check("|v2out(r=300) - v2(long)| / v0", RelDiff(a300.GetV2Out(), f2, v0), 1e-10);
  printf("  raw at r = 300 fm : |v1 - v1(long)| / v0 = %.3e, |v2 - v2(long)| / v0 = %.3e\n",
         RelDiff(raw1, f1, v0), RelDiff(raw2, f2, v0));
  Check("momentum, extrapolated out vs initial (relative)",
        RelDiff(a0.GetV1Out()*m1 + a0.GetV2Out()*m2, P0, P0.Mag()), 1e-14);
  Check("momentum, extrapolated in vs initial (relative)",
        RelDiff(a0.GetV1In()*m1 + a0.GetV2In()*m2, P0, P0.Mag()), 1e-14);
  Check("momentum, long integration vs initial (relative)", RelDiff(Pf, P0, P0.Mag()), 1e-12);
  Check("energy, extrapolated kinetic out vs initial KE+PE (relative)",
        RelDiff(0.5*m1*a0.GetV1Out().Mag2() + 0.5*m2*a0.GetV2Out().Mag2(), E0), 1e-12);
  Check("Vcm at t=0 vs Vcm at r=300 (relative)", RelDiff(a0.GetVcm(), a300.GetVcm(), a0.GetVcm().Mag()), 1e-12);
  if(angles){
    Check("thetaLab1 vs LabAngle(thetaCM, m1, m2) (rad)",
          TMath::Abs(a0.GetThetaLab1() - asymptote::LabAngle(a0.GetThetaCM(), m1, m2)), 1e-12);
    Check("thetaLab2 vs RecoilAngle(thetaCM) (rad)",
          TMath::Abs(a0.GetThetaLab2() - asymptote::RecoilAngle(a0.GetThetaCM())), 1e-12);
    Check("thetaLab1 vs long integration (rad)",
          TMath::Abs(a0.GetThetaLab1() - asymptote::Angle(a0.GetV1In(), f1)), 1e-10);
    if(m1==m2){
      Check("equal masses : |thetaLab1 + thetaLab2 - 90 deg| (rad)",
            TMath::Abs(a0.GetThetaLab1() + a0.GetThetaLab2() - 0.5*TMath::Pi()), 1e-12);
    }
  }else if(m1==m2){
    Check("equal masses head-on : |v1 out| / v0 (projectile stops)", a0.GetV1Out().Mag()/v0, 1e-12);
    Check("equal masses head-on : |v2 out - v1 in| / v0", RelDiff(a0.GetV2Out(), a0.GetV1In(), v0), 1e-12);
  }
}

////////////////////////////////////////////
// 5. head-on details and special states
////////////////////////////////////////////
Bool_t AllFinite(const asymptote & a){
  Double_t vals[] = {a.GetEnergy(), a.GetLMag(), a.GetLRL().Mag(), a.GetEccentricity(),
                     a.GetClosestApproach(), a.GetVinf(), a.GetVInRel().Mag(), a.GetVOutRel().Mag(),
                     a.GetThetaCM(), a.GetImpactParameter(), a.GetV1Out().Mag(), a.GetV2Out().Mag(),
                     a.GetThetaLab1(), a.GetThetaLab2(), a.GetPeriapsisDirection().Mag()};
  for(size_t i=0; i<sizeof(vals)/sizeof(vals[0]); i++) if(!TMath::Finite(vals[i])) return kFALSE;
  return kTRUE;
}

void Special(){
  const Double_t kappa = kCoulombK*2.*79.;
  const Double_t mu = 4.;
  std::cout<<std::endl<<"== 5. Head-on details and special states (alpha, fixed Au) =="<<std::endl;

  asymptote h0(TVector3(-1000., 0., 0.), TVector3(0.05, 0., 0.), mu, kappa);
  Check("head-on : |thetaCM - 180 deg| (rad)", TMath::Abs(h0.GetThetaCM() - TMath::Pi()), 0.);
  Check("head-on : rmin vs HeadOnClosestApproach(v_inf) (relative)",
        RelDiff(h0.GetClosestApproach(), asymptote::HeadOnClosestApproach(h0.GetVinf(), mu, kappa)), 1e-13);
  Check("head-on : rmin vs kappa / E (relative)", RelDiff(h0.GetClosestApproach(), kappa/h0.GetEnergy()), 1e-15);
  Info("head-on : nominal 2 kappa / (mu v0^2) (fm)", asymptote::HeadOnClosestApproach(0.05, mu, kappa));
  Info("head-on : rmin from the asymptote (fm)", h0.GetClosestApproach());
  Check("head-on : RutherfordAngle(b = 0) - 180 deg (rad)",
        TMath::Abs(asymptote::RutherfordAngle(0., 0.05, mu, kappa) - TMath::Pi()), 1e-15);

  // continuity toward L = 0
  const Double_t bs[] = {1e-3, 1e-6, 1e-9};
  for(Int_t i=0; i<3; i++){
    asymptote hb(TVector3(-1000., bs[i], 0.), TVector3(0.05, 0., 0.), mu, kappa);
    Double_t thR = asymptote::RutherfordAngle(hb.GetImpactParameter(), hb.GetVinf(), mu, kappa);
    printf("  b = %6.0e fm : %-8s thetaCM = %.12f deg, rmin = %.9f fm, vout = (%.6e, %.3e, 0)\n",
           bs[i], hb.GetStatusName().Data(), hb.GetThetaCM()*TMath::RadToDeg(), hb.GetClosestApproach(),
           hb.GetVOutRel().X(), hb.GetVOutRel().Y());
    TString name = Form("near head-on b = %.0e : |thetaCM - Rutherford| (rad)", bs[i]);
    Check(name.Data(), TMath::Abs(hb.GetThetaCM() - thR), 1e-12);
    name = Form("near head-on b = %.0e : rmin vs head-on (relative)", bs[i]);
    Check(name.Data(), RelDiff(hb.GetClosestApproach(), h0.GetClosestApproach()), 1e-6);
  }

  // special states
  struct cas{ const char* name; TVector3 r, v; Double_t mu, kappa; asymptote::status expect; };
  const Double_t vpar = TMath::Sqrt(2.*kappa/(mu*1000.));
  cas cs[] = {
    {"repulsive at rest (head-on)",       TVector3(1000.,0.,0.), TVector3(0.,0.,0.),      mu,  kappa, asymptote::kHeadOn},
    {"attractive, bound",                 TVector3(1000.,0.,0.), TVector3(0.,0.001,0.),   mu, -kappa, asymptote::kBound},
    {"attractive, parabolic",             TVector3(1000.,0.,0.), TVector3(0.,vpar,0.),    mu, -kappa, asymptote::kParabolic},
    {"attractive, head-on (radial)",      TVector3(1000.,0.,0.), TVector3(-0.05,0.,0.),   mu, -kappa, asymptote::kRadialAttractive},
    {"attractive, head-on, bound",        TVector3(1000.,0.,0.), TVector3(-0.001,0.,0.),  mu, -kappa, asymptote::kRadialAttractive},
    {"no interaction (kappa = 0)",        TVector3(-1000.,30.,0.), TVector3(0.05,0.,0.),  mu,  0.,    asymptote::kFree},
    {"zero distance",                     TVector3(0.,0.,0.),    TVector3(0.05,0.,0.),    mu,  kappa, asymptote::kInvalid},
    {"zero mass",                         TVector3(1000.,0.,0.), TVector3(0.05,0.,0.),    0.,  kappa, asymptote::kInvalid},
  };
  for(size_t i=0; i<sizeof(cs)/sizeof(cs[0]); i++){
    asymptote a(cs[i].r, cs[i].v, cs[i].mu, cs[i].kappa);
    TString name = Form("%-30s -> %-17s E %10.3e rmin %9.3f", cs[i].name, a.GetStatusName().Data(),
                        a.GetEnergy(), a.GetClosestApproach());
    CheckTrue(name.Data(), a.GetStatus()==cs[i].expect && AllFinite(a));
  }
  {
    asymptote a(TVector3(-1000.,30.,0.), TVector3(0.05,0.,0.), mu, 0.);
    Check("kappa = 0 : thetaCM (rad)", a.GetThetaCM(), 0.);
    Check("kappa = 0 : |b - 30| and |rmin - 30| (fm)",
          TMath::Max(TMath::Abs(a.GetImpactParameter()-30.), TMath::Abs(a.GetClosestApproach()-30.)), 1e-12);
  }
  {
    asymptote a(TVector3(-1000.,0.,0.), TVector3(0.05,0.,0.), 4., kTRUE,
                TVector3(), TVector3(), 206., kTRUE, kappa);
    CheckTrue("both fixed -> INVALID", a.GetStatus()==asymptote::kInvalid && AllFinite(a));
  }
  {
    // the order of the particles does not matter (fixed centre given first)
    asymptote a12(TVector3(-1000.,50.,0.), TVector3(0.05,0.,0.), 4., kFALSE,
                  TVector3(), TVector3(), 206., kTRUE, kappa);
    asymptote a21(TVector3(), TVector3(), 206., kTRUE,
                  TVector3(-1000.,50.,0.), TVector3(0.05,0.,0.), 4., kFALSE, kappa);
    Check("swapped order : |v1out(12) - v2out(21)| / v0",
          RelDiff(a12.GetV1Out(), a21.GetV2Out(), 0.05), 1e-15);
    Check("swapped order : |thetaLab1(12) - thetaLab2(21)| (rad)",
          TMath::Abs(a12.GetThetaLab1() - a21.GetThetaLab2()), 1e-15);
  }
}

int main(int argc, char** argv){
  Invariance();
  FixedTarget("2-3. Repulsive, fixed target", 4., 2., 206., 79., 10.,  0.05);
  FixedTarget("2-3. Repulsive, fixed target", 4., 2., 206., 79., 50.,  0.05);
  FixedTarget("2-3. Repulsive, fixed target", 4., 2., 206., 79., 200., 0.05);
  FixedTarget("2-3. Head-on, fixed target",   4., 2., 206., 79., 0.,   0.05);
  FixedTarget("2-3. Attractive (opposite charges), fixed target", 4., -2., 206., 79., 50., 0.05);
  FreeTarget("4. Free target, alpha on alpha", 4., 2., 4., 2., 2., 0.05, kTRUE);
  FreeTarget("4. Free target, alpha on Au",    4., 2., 206., 79., 20., 0.05, kTRUE);
  FreeTarget("4. Free target, alpha on alpha, head-on", 4., 2., 4., 2., 0., 0.05, kFALSE);
  Special();

  std::cout<<std::endl<<"== "<<(gCheck-gFail)<<" / "<<gCheck<<" checks passed =="<<std::endl;
  return gFail==0 ? 0 : 1;
}
