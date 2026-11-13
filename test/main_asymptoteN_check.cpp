////////////////////////////////////////////
//
// RooParticle Test Application
// N-body asymptotic velocities (asymptoteN.hpp)
//
// 1. two bodies, pure Coulomb : every method against the exact class
//    asymptote at several stopping distances ; measured order in 1/r.
// 2. three and four bodies : the corrected asymptotes at several stopping
//    distances against the brute-force reference (computed from the same
//    state) ; actual error / reported estimate ; cost.
// 3. the residual (remaining deflection) as a termination criterion :
//    monotonic once separating, and its ratio to the true remaining deflection.
// 4. incoming asymptotes of a run started at a finite distance.
// 5. energy and momentum at infinity.
// 6. bound subsystems : a bound pair, a particle bound to a fixed centre.
// 7. Coulomb + Woods-Saxon (short range) : two and three bodies.
// 8. relativistic mode : first order against the relativistic brute force.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "asymptote.hpp"
#include "asymptoteN.hpp"

#include "iostream"
#include "vector"
#include <cstdio>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"

static const Double_t kCoulombK = 0.00154629;

static coulombForce    gCoulomb;
static woodsSaxonForce gWS;       // defaults : 50 MeV, R = 8.886 fm, a = 0.65 fm (alpha + Au)
static Int_t gFail = 0;
static Int_t gCheck = 0;
static Double_t gWorstRatio = 0.; // largest actual / estimated error, two bodies (sections 1, 7)
static Double_t gWorstPairN = 0.; // pairwise and auto, three and four bodies (section 2)
static Double_t gWorstBornN = 0.; // Born, three and four bodies (section 2)

void Check(const char* what, Double_t val, Double_t tol){
  Bool_t ok = TMath::Finite(val) && val <= tol;
  printf("  %-66s %11.3e < %8.1e  %s\n", what, val, tol, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
void CheckRange(const char* what, Double_t val, Double_t lo, Double_t hi){
  Bool_t ok = TMath::Finite(val) && val >= lo && val <= hi;
  printf("  %-66s %11.4f in [%.2f, %.2f]  %s\n", what, val, lo, hi, ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
void CheckTrue(const char* what, Bool_t ok){
  printf("  %-66s %24s  %s\n", what, "", ok ? "ok" : "FAIL");
  gCheck++;
  if(!ok) gFail++;
}
void Info(const char* what, Double_t val, const char* unit=""){
  printf("  %-66s %16.8g %s\n", what, val, unit);
}

////////////////////////////////////////////
// helpers
////////////////////////////////////////////
event* NewEvent(Bool_t nuclear = kFALSE){
  event* ev = new event();
  ev->AddForce(&gCoulomb);
  if(nuclear) ev->AddForce(&gWS);
  ev->SetIntegrator(event::kRK45);
  ev->SetTolerance(1e-13, 1e-15);
  return ev;
}

// every pair (not both fixed) receding and farther than D
Bool_t Separated(event* ev, Double_t D){
  const Int_t n = ev->getNParticle();
  for(Int_t i=0; i<n; i++){
    for(Int_t j=i+1; j<n; j++){
      particle* a = ev->getParticle(i);
      particle* b = ev->getParticle(j);
      if(a->IsInvincible() && b->IsInvincible()) continue;
      TVector3 r = a->GetX() - b->GetX();
      TVector3 w = (a->IsInvincible() ? TVector3() : a->GetV()) - (b->IsInvincible() ? TVector3() : b->GetV());
      if(r.Dot(w) < 0. || r.Mag() < D) return kFALSE;
    }
  }
  return kTRUE;
}

// Integrate (RK45, tolerance of the event) until Separated(D). Chunks are
// small compared to the distances, so the stop is close to D.
void IntegrateTo(event* ev, Double_t D){
  for(Long_t guard=0; guard<10000000; guard++){
    if(Separated(ev, D)) return;
    Double_t dmin = ev->GetMinDistance();
    Double_t vmax = 1e-6;
    for(Int_t i=0; i<ev->getNParticle(); i++) vmax = TMath::Max(vmax, ev->getParticle(i)->GetV().Mag());
    ev->DeriveDT(TMath::Max(0.2, 0.02*dmin/(2.*vmax)));
  }
  std::cout<<"  IntegrateTo did not end"<<std::endl;
  gFail++;
}

Double_t MinPairDistance(event* ev){
  Double_t d = 1e300;
  for(Int_t i=0; i<ev->getNParticle(); i++){
    for(Int_t j=i+1; j<ev->getNParticle(); j++){
      if(ev->getParticle(i)->IsInvincible() && ev->getParticle(j)->IsInvincible()) continue;
      d = TMath::Min(d, (ev->getParticle(i)->GetX() - ev->getParticle(j)->GetX()).Mag());
    }
  }
  return d;
}

// largest |v_inf(a) - v_ref| over the moving particles
Double_t ErrVs(const asymptoteN & a, const std::vector<TVector3> & vref){
  Double_t e = 0.;
  for(Int_t i=0; i<a.GetNParticle(); i++){
    if(a.IsFixed(i)) continue;
    e = TMath::Max(e, (a.GetVInf(i) - vref[i]).Mag());
  }
  return e;
}
// largest actual / (estimate + reference error) over the moving particles
Double_t RatioVs(const asymptoteN & a, const std::vector<TVector3> & vref, const std::vector<Double_t> & eref){
  Double_t r = 0.;
  for(Int_t i=0; i<a.GetNParticle(); i++){
    if(a.IsFixed(i)) continue;
    r = TMath::Max(r, (a.GetVInf(i) - vref[i]).Mag()/(a.GetError(i) + eref[i]));
  }
  return r;
}
std::vector<TVector3> VInf(const asymptoteN & a){
  std::vector<TVector3> v;
  for(Int_t i=0; i<a.GetNParticle(); i++) v.push_back(a.GetVInf(i));
  return v;
}
std::vector<Double_t> Errs(const asymptoteN & a){
  std::vector<Double_t> v;
  for(Int_t i=0; i<a.GetNParticle(); i++) v.push_back(a.GetError(i));
  return v;
}
Double_t Order(Double_t e1, Double_t e2, Double_t d1, Double_t d2){
  return TMath::Log(e1/e2)/TMath::Log(d2/d1);
}
Double_t KineticSum(event* ev, const asymptoteN & a){
  Double_t k = 0.;
  for(Int_t i=0; i<a.GetNParticle(); i++){
    if(a.IsFixed(i)) continue;
    k += 0.5*ev->getParticle(i)->GetM0()*a.GetVInf(i).Mag2();
  }
  return k;
}
Double_t NetE(event* ev){ return ev->GetNetKE() + ev->GetNetPE(); } // at most one fixed particle here
// energy uncertainty implied by the velocity error estimates : sum m |v_inf| err
Double_t EnergyAllowance(event* ev, const asymptoteN & a){
  Double_t s = 0.;
  for(Int_t i=0; i<a.GetNParticle(); i++){
    if(a.IsFixed(i)) continue;
    s += ev->getParticle(i)->GetM0()*a.GetVInf(i).Mag()*a.GetError(i);
  }
  return s;
}

////////////////////////////////////////////
// 1. two bodies, pure Coulomb, against the exact asymptote
////////////////////////////////////////////
void TwoBody(const char* title, Bool_t fixedTarget, Double_t m2, Double_t q2, Double_t b0){
  const Double_t m1 = 4., q1 = 2., v0 = 0.05;
  const Double_t kappa = kCoulombK*q1*q2;
  printf("\n== 1. %s (alpha on %s m %.0f q %.0f, b = %.0f fm, v = %.2f c, from 1000 fm) ==\n",
         title, fixedTarget ? "fixed" : "free", m2, q2, b0, v0);
  event* ev = NewEvent();
  ev->AddParticle(new EMparticle(m1, q1, TVector3(-1000., b0, 0.), TVector3(v0, 0., 0.)));
  ev->AddParticle(new EMparticle(m2, q2, TVector3(), TVector3(), fixedTarget));

  const Int_t nd = 6;
  const Double_t Ds[nd] = {200., 500., 1000., 3000., 10000., 30000.};
  Double_t eRaw[nd], eBorn[nd], ePair[nd], eGen[nd], eGenB[nd], eAuto[nd], eBrute[nd];
  Double_t worst = 0., worstBrute = 0., worstAuto = 0.;
  printf("  %7s %9s | %9s | %9s %9s | %9s %9s | %9s %9s %9s | %9s %9s | %9s %9s %6s\n",
         "D(fm)", "|PE/KE|", "raw", "Born", "est", "pair", "est", "gen.pair", "est", "gen.Born",
         "auto", "est", "brute", "est", "subst");
  for(Int_t k=0; k<nd; k++){
    IntegrateTo(ev, Ds[k]);
    particle* p1 = ev->getParticle(0);
    particle* p2 = ev->getParticle(1);
    asymptote ex(p1, p2, kappa);
    std::vector<TVector3> vex(2);
    vex[0] = ex.GetV1Out();
    vex[1] = fixedTarget ? TVector3() : ex.GetV2Out();
    std::vector<Double_t> zero(2, 0.);

    asymptoteN born(ev, asymptoteN::kBorn);
    asymptoteN pair(ev, asymptoteN::kPairwise);
    asymptoteN gen;  gen.SetAnalytic(kFALSE);  gen.Compute(ev, asymptoteN::kPairwise);
    asymptoteN genb; genb.SetAnalytic(kFALSE); genb.Compute(ev, asymptoteN::kBorn);
    asymptoteN aut(ev, asymptoteN::kAuto);
    asymptoteN bru(ev, asymptoteN::kBrute);

    Double_t raw = 0.;
    raw = TMath::Max((p1->GetV() - vex[0]).Mag(), fixedTarget ? 0. : (p2->GetV() - vex[1]).Mag());
    eRaw[k]   = raw/v0;
    eBorn[k]  = ErrVs(born, vex)/v0;
    ePair[k]  = ErrVs(pair, vex)/v0;
    eGen[k]   = ErrVs(gen,  vex)/v0;
    eGenB[k]  = ErrVs(genb, vex)/v0;
    eAuto[k]  = ErrVs(aut,  vex)/v0;
    eBrute[k] = ErrVs(bru,  vex)/v0;
    Double_t r = 0.;
    r = TMath::Max(r, RatioVs(born, vex, zero));
    r = TMath::Max(r, RatioVs(pair, vex, zero));
    r = TMath::Max(r, RatioVs(gen,  vex, zero));
    r = TMath::Max(r, RatioVs(genb, vex, zero));
    worst      = TMath::Max(worst, r);
    worstAuto  = TMath::Max(worstAuto, RatioVs(aut, vex, zero));
    worstBrute = TMath::Max(worstBrute, RatioVs(bru, vex, zero));
    Double_t ke = 0.5*m1*p1->GetV().Mag2() + (fixedTarget ? 0. : 0.5*m2*p2->GetV().Mag2());
    printf("  %7.0f %9.2e | %9.2e | %9.2e %9.2e | %9.2e %9.2e | %9.2e %9.2e %9.2e | %9.2e %9.2e | %9.2e %9.2e %6ld\n",
           Ds[k], TMath::Abs(ev->GetNetPE())/ke, eRaw[k], eBorn[k], born.GetMaxError()/v0,
           ePair[k], pair.GetMaxError()/v0, eGen[k], gen.GetMaxError()/v0, eGenB[k],
           eAuto[k], aut.GetMaxError()/v0, eBrute[k], bru.GetMaxError()/v0, bru.GetNSubsteps());
  }
  printf("  (errors : max over the moving particles of |v_inf - exact| / v0 ; est : reported estimate / v0)\n");
  printf("  measured order in 1/D :  %7s %7s %7s\n", "raw", "Born", "genBorn");
  for(Int_t k=0; k+1<nd; k++){
    printf("    D %6.0f -> %6.0f :  %7.3f %7.3f %7.3f\n", Ds[k], Ds[k+1],
           Order(eRaw[k], eRaw[k+1], Ds[k], Ds[k+1]), Order(eBorn[k], eBorn[k+1], Ds[k], Ds[k+1]),
           Order(eGenB[k], eGenB[k+1], Ds[k], Ds[k+1]));
  }
  Double_t mPair = 0., mGen = 0., mAuto = 0., mBrute = 0.;
  for(Int_t k=0; k<nd; k++){
    mPair = TMath::Max(mPair, ePair[k]); mGen = TMath::Max(mGen, eGen[k]);
    mAuto = TMath::Max(mAuto, eAuto[k]); mBrute = TMath::Max(mBrute, eBrute[k]);
  }
  CheckRange("order of the raw velocity error (last interval)", Order(eRaw[nd-2], eRaw[nd-1], Ds[nd-2], Ds[nd-1]), 0.9, 1.1);
  CheckRange("order of the first-order (Born) error (last interval)", Order(eBorn[nd-2], eBorn[nd-1], Ds[nd-2], Ds[nd-1]), 1.8, 2.2);
  CheckRange("order of the quadrature Born error (last interval)", Order(eGenB[nd-2], eGenB[nd-1], Ds[nd-2], Ds[nd-1]), 1.8, 2.2);
  Check("pairwise (closed form) : max error / v0 (exact for two bodies)", mPair, 1e-14);
  Check("pairwise (orbit quadrature) : max error / v0", mGen, 1e-13);
  Check("auto : max error / v0", mAuto, 1e-13);
  Check("brute force : max error / v0", mBrute, 1e-11);
  Check("Born, pairwise : max actual / estimated error", worst, 1.);
  Check("auto : max actual / estimated error", worstAuto, 1.);
  Check("brute force : max actual / estimated error", worstBrute, 1.);
  gWorstRatio = TMath::Max(gWorstRatio, worst);
  delete ev;
}

////////////////////////////////////////////
// 2. N bodies against the brute force
////////////////////////////////////////////
struct spec { Double_t m, q; TVector3 x, v; Bool_t fixed; };

event* Build(const std::vector<spec> & s, Bool_t nuclear = kFALSE){
  event* ev = NewEvent(nuclear);
  for(size_t i=0; i<s.size(); i++) ev->AddParticle(new EMparticle(s[i].m, s[i].q, s[i].x, s[i].v, s[i].fixed));
  return ev;
}

// fixed Au + two alphas arriving together from perpendicular directions
std::vector<spec> CaseFixed3(){
  std::vector<spec> s;
  s.push_back({197., 79., TVector3(), TVector3(), kTRUE});
  s.push_back({4., 2., TVector3(-1000., 40., 0.), TVector3(0.05, 0., 0.), kFALSE});
  s.push_back({4., 2., TVector3(30., -1000., 10.), TVector3(0., 0.05, 0.), kFALSE});
  return s;
}
// three free bodies : C12 at rest, two alphas
std::vector<spec> CaseFree3(){
  std::vector<spec> s;
  s.push_back({12., 6., TVector3(), TVector3(), kFALSE});
  s.push_back({4., 2., TVector3(-1000., 15., 0.), TVector3(0.05, 0., 0.), kFALSE});
  s.push_back({4., 2., TVector3(10., -1000., -8.), TVector3(0., 0.05, 0.), kFALSE});
  return s;
}
// four free bodies : Au at rest, two alphas and a proton
std::vector<spec> CaseFree4(){
  std::vector<spec> s;
  s.push_back({197., 79., TVector3(), TVector3(), kFALSE});
  s.push_back({4., 2., TVector3(-1000., 30., 0.), TVector3(0.05, 0., 0.), kFALSE});
  s.push_back({4., 2., TVector3(20., -1000., 15.), TVector3(0., 0.05, 0.), kFALSE});
  s.push_back({1., 1., TVector3(-40., 25., -1000.), TVector3(0., 0., 0.05), kFALSE});
  return s;
}

void NBody(const char* title, const std::vector<spec> & s, Bool_t nuclear = kFALSE,
           Double_t dFirst = 200., Double_t v0 = 0.05){
  printf("\n== 2. %s ==\n", title);
  event* ev = Build(s, nuclear);
  Bool_t anyFixed = kFALSE;
  for(size_t i=0; i<s.size(); i++) anyFixed |= s[i].fixed;
  const Double_t E0 = NetE(ev);
  const TVector3 P0 = ev->GetNetMomentum();

  std::vector<Double_t> Ds;
  for(Double_t d = dFirst; d <= 40000.; d *= 2.5) Ds.push_back(d);
  const Int_t nd = Ds.size();
  std::vector<Double_t> eBorn(nd), ePair(nd), eAuto(nd), estB(nd), estP(nd);
  std::vector<TVector3> vFirst;
  Double_t worstB = 0., worstP = 0., worstA = 0., refSpread = 0., eDefB = 0., eDefP = 0.;
  Long_t subMax = 0;
  Double_t tBorn = 0., tPair = 0., tAuto = 0., tBrute = 0.;
  printf("  %7s %8s %9s | %9s %9s %6s | %9s %9s %6s | %9s %9s %3s | %9s %6s %s\n",
         "D(fm)", "eta_max", "raw", "Born", "est", "ratio", "pair", "est", "ratio", "auto", "est", "blk",
         "brute est", "subst", "flags(pair)");
  for(Int_t k=0; k<nd; k++){
    IntegrateTo(ev, Ds[k]);
    asymptoteN bru(ev, asymptoteN::kBrute);
    std::vector<TVector3> vref = VInf(bru);
    std::vector<Double_t> eref = Errs(bru);
    if(k==0) vFirst = vref;
    refSpread = TMath::Max(refSpread, ErrVs(bru, vFirst)/v0);
    asymptoteN born(ev, asymptoteN::kBorn);
    asymptoteN pair(ev, asymptoteN::kPairwise);
    asymptoteN aut(ev, asymptoteN::kAuto);
    Double_t raw = 0.;
    for(Int_t i=0; i<ev->getNParticle(); i++){
      if(s[i].fixed) continue;
      raw = TMath::Max(raw, (ev->getParticle(i)->GetV() - vref[i]).Mag());
    }
    eBorn[k] = ErrVs(born, vref)/v0;  estB[k] = born.GetMaxError()/v0;
    ePair[k] = ErrVs(pair, vref)/v0;  estP[k] = pair.GetMaxError()/v0;
    eAuto[k] = ErrVs(aut, vref)/v0;
    Double_t rB = RatioVs(born, vref, eref), rP = RatioVs(pair, vref, eref), rA = RatioVs(aut, vref, eref);
    worstB = TMath::Max(worstB, rB);
    worstP = TMath::Max(worstP, rP);
    worstA = TMath::Max(worstA, rA);
    subMax = TMath::Max(subMax, bru.GetNSubsteps());
    tBorn += born.GetRunTime(); tPair += pair.GetRunTime(); tAuto += aut.GetRunTime(); tBrute += bru.GetRunTime();
    eDefB = TMath::Max(eDefB, TMath::Abs(born.GetEnergyInf() - born.GetEnergyNow())/born.GetEnergyNow());
    eDefP = TMath::Max(eDefP, TMath::Abs(pair.GetEnergyInf() - pair.GetEnergyNow())/pair.GetEnergyNow());
    printf("  %7.0f %8.2e %9.2e | %9.2e %9.2e %6.3f | %9.2e %9.2e %6.3f | %9.2e %9.2e %3d | %9.2e %6ld %s\n",
           Ds[k], pair.GetMaxEta(), raw/v0, eBorn[k], estB[k], rB, ePair[k], estP[k], rP,
           eAuto[k], aut.GetMaxError()/v0, aut.GetNBlocks(), bru.GetMaxError()/v0, bru.GetNSubsteps(),
           pair.GetFlagNames(pair.GetFlags()).Data());
    if(k==nd-1){
      // 5. energy and momentum at infinity (brute force and pairwise, last stop)
      Check("brute force : |sum KE_inf - E_total| / E_total", TMath::Abs(KineticSum(ev, bru) - E0)/E0, 1e-10);
      Check("pairwise    : |sum KE_inf - E_total| / (sum m |v_inf| err)",
            TMath::Abs(KineticSum(ev, pair) - E0)/EnergyAllowance(ev, pair), 1.);
      Info("pairwise    : |sum KE_inf - E_total| / E_total", TMath::Abs(KineticSum(ev, pair) - E0)/E0);
      if(!anyFixed){
        Check("brute force : |P_inf - P_0| / |P_0|", (bru.GetMomentumInf() - P0).Mag()/P0.Mag(), 1e-10);
        Check("pairwise    : |P_inf - P_0| / |P_0| (exact by construction)", (pair.GetMomentumInf() - P0).Mag()/P0.Mag(), 1e-13);
      }
    }
  }
  printf("  (errors / v0 against the brute force from the same state ; ratio = actual / (estimate + reference estimate))\n");
  printf("  measured order in 1/D (pairs of stops) :");
  for(Int_t k=0; k+1<nd; k++) printf(" B %.2f P %.2f |", Order(eBorn[k], eBorn[k+1], Ds[k], Ds[k+1]), Order(ePair[k], ePair[k+1], Ds[k], Ds[k+1]));
  printf("\n");
  Info("cost per evaluation (ms) : Born", 1e3*tBorn/nd);
  Info("cost per evaluation (ms) : pairwise", 1e3*tPair/nd);
  Info("cost per evaluation (ms) : auto", 1e3*tAuto/nd);
  Info("cost per evaluation (ms) : brute force", 1e3*tBrute/nd);
  Info("largest number of RK45 substeps of the brute force", (Double_t) subMax);
  Info("spread of the brute-force reference along the trajectory / v0", refSpread);
  Info("energy defect of Born at infinity (max, relative)", eDefB);
  Info("energy defect of pairwise at infinity (max, relative)", eDefP);
  Check("pairwise : error at the last stop / error at the first stop", ePair[nd-1]/ePair[0], 1e-3);
  Check("Born     : error at the last stop / error at the first stop", eBorn[nd-1]/eBorn[0], 1e-3);
  Check("Born     : max actual / estimated error", worstB, 1.);
  Check("pairwise : max actual / estimated error", worstP, 1.);
  Check("auto     : max actual / estimated error", worstA, 1.);
  Check("brute-force reference consistent along the trajectory (/ v0)", refSpread, 1e-9);
  gWorstPairN = TMath::Max(gWorstPairN, TMath::Max(worstP, worstA));
  gWorstBornN = TMath::Max(gWorstBornN, worstB);
  delete ev;
}

////////////////////////////////////////////
// 3. residual as termination criterion
////////////////////////////////////////////
void Residual(const char* title, const std::vector<spec> & s){
  printf("\n== 3. Residual (remaining deflection) as a stop criterion : %s ==\n", title);
  event* ev = Build(s);
  IntegrateTo(ev, 0.);                              // every pair receding
  asymptoteN bru(ev, asymptoteN::kBrute);           // the trajectory's true v_inf
  std::vector<TVector3> vref = VInf(bru);
  Info("closest pair distance when all pairs start to recede (fm)", MinPairDistance(ev));

  Double_t prev = 1e300, worstIncrease = 0.;
  Double_t rbMin = 1e300, rbMax = 0., rpMin = 1e300, rpMax = 0.;
  Int_t nsample = 0;
  Bool_t allSep = kTRUE;
  printf("  %8s %11s %11s %8s %8s %8s %8s\n", "D(fm)", "predicted", "actual", "Born/act", "", "pair/act", "");
  printf("  %8s %11s %11s %8s %8s %8s %8s\n", "", "max (rad)", "max (rad)", "min", "max", "min", "max");
  Double_t nextPrint = 100.;
  while(MinPairDistance(ev) < 30000.){
    Bool_t sep;
    const Double_t pred = asymptoteN::RemainingDeflection(ev, &sep);
    allSep &= sep;
    if(pred > prev) worstIncrease = TMath::Max(worstIncrease, (pred - prev)/prev);
    prev = pred;
    nsample++;
    const Double_t d = MinPairDistance(ev);
    if(d >= nextPrint){
      asymptoteN born; born.SetErrorEstimate(kFALSE); born.Compute(ev, asymptoteN::kBorn);
      asymptoteN pair(ev, asymptoteN::kPairwise);
      Double_t act = 0., bmin = 1e300, bmax = 0., pmin = 1e300, pmax = 0.;
      for(Int_t i=0; i<ev->getNParticle(); i++){
        if(s[i].fixed) continue;
        const Double_t a = asymptote::Angle(ev->getParticle(i)->GetV(), vref[i]);
        act = TMath::Max(act, a);
        bmin = TMath::Min(bmin, born.GetDeflection(i)/a); bmax = TMath::Max(bmax, born.GetDeflection(i)/a);
        pmin = TMath::Min(pmin, pair.GetDeflection(i)/a); pmax = TMath::Max(pmax, pair.GetDeflection(i)/a);
      }
      printf("  %8.0f %11.4e %11.4e %8.4f %8.4f %8.5f %8.5f\n", d, pred, act, bmin, bmax, pmin, pmax);
      rbMin = TMath::Min(rbMin, bmin); rbMax = TMath::Max(rbMax, bmax);
      rpMin = TMath::Min(rpMin, pmin); rpMax = TMath::Max(rpMax, pmax);
      nextPrint *= 2.;
    }
    // about 40 samples per doubling of the distance
    Double_t vmax = 0.;
    for(Int_t i=0; i<ev->getNParticle(); i++) vmax = TMath::Max(vmax, ev->getParticle(i)->GetV().Mag());
    ev->DeriveDT(0.0175*d/vmax);
  }
  printf("  (%d samples ; predicted = RemainingDeflection(ev), first order ; actual = angle(v, v_inf brute))\n", nsample);
  CheckTrue("every interacting pair reported receding (IsSeparating)", allSep);
  Check("largest relative increase between successive samples", worstIncrease, 1e-12);
  CheckRange("first-order prediction / actual remaining deflection : min", rbMin, 0.5, 2.);
  CheckRange("first-order prediction / actual remaining deflection : max", rbMax, 0.5, 2.);
  CheckRange("pairwise prediction / actual remaining deflection : min", rpMin, 0.8, 1.25);
  CheckRange("pairwise prediction / actual remaining deflection : max", rpMax, 0.8, 1.25);
  delete ev;
}

////////////////////////////////////////////
// 4. incoming asymptotes
////////////////////////////////////////////
void Incoming(const char* title, const std::vector<spec> & s){
  printf("\n== 4. Incoming asymptotes of a run started at 1000 fm : %s ==\n", title);
  event* ev = Build(s);
  Bool_t anyFixed = kFALSE;
  for(size_t i=0; i<s.size(); i++) anyFixed |= s[i].fixed;
  const Double_t E0 = NetE(ev);
  const TVector3 P0 = ev->GetNetMomentum();
  asymptoteN inP(ev, asymptoteN::kPairwise, asymptoteN::kIncoming);
  asymptoteN inB(ev, asymptoteN::kBorn,     asymptoteN::kIncoming);
  asymptoteN inA(ev, asymptoteN::kAuto,     asymptoteN::kIncoming);
  asymptoteN inR(ev, asymptoteN::kBrute,    asymptoteN::kIncoming);
  std::vector<TVector3> vref = VInf(inR);
  std::vector<Double_t> eref = Errs(inR);
  printf("  %4s %10s %12s %12s %12s | %10s %10s %10s\n", "part", "|v| nomin.", "|v_in|", "tilt(deg)", "|dv|/v0",
         "pair err", "pair est", "auto err");
  for(Int_t i=0; i<ev->getNParticle(); i++){
    if(s[i].fixed) continue;
    TVector3 v = ev->getParticle(i)->GetV();
    printf("  %4d %10.6f %12.9f %12.6f %12.4e | %10.2e %10.2e %10.2e\n", i, v.Mag(), vref[i].Mag(),
           (v.Mag2()>0. ? asymptote::Angle(v, vref[i])*TMath::RadToDeg() : 0.), (vref[i]-v).Mag()/0.05,
           (inP.GetVInf(i)-vref[i]).Mag(), inP.GetError(i), (inA.GetVInf(i)-vref[i]).Mag());
  }
  Check("pairwise incoming : max actual / estimated error", RatioVs(inP, vref, eref), 1.);
  Check("Born incoming : max actual / estimated error", RatioVs(inB, vref, eref), 1.);
  Check("auto incoming : max actual / estimated error", RatioVs(inA, vref, eref), 1.);
  Check("brute incoming : |sum KE_in - E_total| / E_total", TMath::Abs(KineticSum(ev, inR) - E0)/E0, 1e-10);
  if(!anyFixed) Check("brute incoming : |P_in - P_0| / |P_0|", (inR.GetMomentumInf() - P0).Mag()/P0.Mag(), 1e-10);
  // the same incoming asymptote from after the collision : integrated back through it
  IntegrateTo(ev, 2000.);
  asymptoteN inR2(ev, asymptoteN::kBrute, asymptoteN::kIncoming);
  Info("RK45 substeps of the backward brute force through the collision", (Double_t) inR2.GetNSubsteps());
  Check("incoming from after the collision vs from the start (brute, / v0)", ErrVs(inR2, vref)/0.05, 1e-9);
  asymptoteN out(ev, asymptoteN::kPairwise);
  CheckTrue("from after the collision : outgoing state is separating", out.IsSeparating());
  CheckTrue("from after the collision : the incoming pairwise state is not (collision ahead)", !asymptoteN(ev, asymptoteN::kPairwise, asymptoteN::kIncoming).IsSeparating());
  delete ev;
}

////////////////////////////////////////////
// 6. bound subsystems
////////////////////////////////////////////
void Bound(){
  printf("\n== 6. Bound subsystems ==\n");
  // (a) a bound pair (+1, -1 ; separation 10 fm) moving at 0.02 c, and an alpha passing at b = 200 fm
  event* ev = NewEvent();
  const TVector3 Vcm(0.02, 0., 0.), w(0.012, 0., 0.);
  ev->AddParticle(new EMparticle(1.,  1., TVector3(0.,  5., 0.), Vcm + w*0.5));
  ev->AddParticle(new EMparticle(1., -1., TVector3(0., -5., 0.), Vcm - w*0.5));
  ev->AddParticle(new EMparticle(4.,  2., TVector3(1000., 200., 0.), TVector3(-0.05, 0., 0.)));
  const Double_t E0 = NetE(ev);
  for(Int_t stage=0; stage<2; stage++){
    if(stage==1){
      // until the alpha has passed and is 1000 fm away, receding
      for(Long_t g=0; g<10000000; g++){
        TVector3 cm = (ev->getParticle(0)->GetX() + ev->getParticle(1)->GetX())*0.5;
        TVector3 r = ev->getParticle(2)->GetX() - cm;
        TVector3 vr = ev->getParticle(2)->GetV() - (ev->getParticle(0)->GetV() + ev->getParticle(1)->GetV())*0.5;
        if(r.Dot(vr) > 0. && r.Mag() > 1000.) break;
        ev->DeriveDT(20.);
      }
    }
    asymptoteN pair(ev, asymptoteN::kPairwise);
    asymptoteN bru; bru.SetBruteTolerance(1e-8); bru.Compute(ev, asymptoteN::kBrute);
    printf("  %s : pairwise flags %s, clusters %d", stage==0 ? "at the start      " : "after the passage ",
           pair.GetFlagNames(pair.GetFlags()).Data(), pair.GetNClusters());
    if(pair.GetNClusters()>0) printf(", E_int %.6e, isolation %.2e", pair.GetClusterEnergy(0), pair.GetClusterIsolation(0));
    printf("\n");
    printf("    raw v (+) %.6f (-) %.6f ; cluster V_inf (%.9e, %.3e) ; brute (%.9e, %.3e), %ld substeps\n",
           ev->getParticle(0)->GetV().Mag(), ev->getParticle(1)->GetV().Mag(),
           pair.GetVInf(0).X(), pair.GetVInf(0).Y(), bru.GetVInf(0).X(), bru.GetVInf(0).Y(), bru.GetNSubsteps());
    TString tag = (stage==0) ? "start : " : "after : ";
    CheckTrue((tag + "one cluster {0, 1} reported, the alpha free").Data(),
              pair.GetNClusters()==1 && pair.GetCluster(0)==0 && pair.GetCluster(1)==0 && pair.GetCluster(2)==-1
              && pair.GetClusterMembers(0).size()==2);
    CheckTrue((tag + "members flagged CLUSTER, isolated, result reliable").Data(),
              (pair.GetFlags(0) & asymptoteN::kInCluster) && (pair.GetFlags(1) & asymptoteN::kInCluster)
              && !(pair.GetFlags() & asymptoteN::kNotIsolated) && pair.IsReliable());
    CheckTrue((tag + "both members get the cluster velocity").Data(),
              (pair.GetVInf(0) - pair.GetVInf(1)).Mag()==0. && (pair.GetVInf(0) - pair.GetClusterVInf(0)).Mag()==0.);
    CheckTrue((tag + "brute force finds the same cluster").Data(),
              bru.GetNClusters()==1 && bru.GetCluster(0)==0 && bru.GetCluster(1)==0 && bru.GetCluster(2)==-1);
    Check((tag + "pairwise vs brute : max actual / estimated error").Data(), RatioVs(pair, VInf(bru), Errs(bru)), 1.);
    Check((tag + "KE_inf + E_int vs E_total (pairwise), / (sum m|v|err + |E0| 1e-12)").Data(),
          TMath::Abs(pair.GetEnergyInf() - E0)/(EnergyAllowance(ev, pair) + 1e-12*TMath::Abs(E0)), 1.);
  }
  delete ev;

  // (b) a light particle bound to a fixed Au, and an alpha passing
  event* ev2 = NewEvent();
  ev2->AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), kTRUE));
  ev2->AddParticle(new EMparticle(1., -1., TVector3(100., 0., 0.), TVector3(0., 0.02, 0.)));
  ev2->AddParticle(new EMparticle(4., 2., TVector3(-1000., 300., 0.), TVector3(0.05, 0., 0.)));
  asymptoteN a2(ev2, asymptoteN::kPairwise);
  printf("  bound to a fixed centre : flags of the light particle %s, of the alpha %s\n",
         a2.GetFlagNames(a2.GetFlags(1)).Data(), a2.GetFlagNames(a2.GetFlags(2)).Data());
  CheckTrue("bound to fixed : flagged BOUND_TO_FIXED, v_inf = 0, error 0",
            (a2.GetFlags(1) & asymptoteN::kBoundToFixed) && a2.GetVInf(1).Mag()==0. && a2.GetError(1)==0.);
  // the alpha sees the composite charge 79 - 1 at the fixed centre
  asymptote ex(ev2->getParticle(2)->GetX(), ev2->getParticle(2)->GetV(), 4., kFALSE,
               TVector3(), TVector3(), 197., kTRUE, kCoulombK*2.*78.);
  Info("alpha : |v_inf - two-body(q = 78)| / v0", (a2.GetVInf(2) - ex.GetV1Out()).Mag()/0.05);
  Info("alpha : error estimate / v0 (multipole of the orbiting charge)", a2.GetError(2)/0.05);
  Check("alpha : |v_inf - two-body(q = 78)| / estimate", (a2.GetVInf(2) - ex.GetV1Out()).Mag()/a2.GetError(2), 1.);
  delete ev2;
}

////////////////////////////////////////////
// 7. Coulomb + Woods-Saxon
////////////////////////////////////////////
void Nuclear2(){
  const Double_t v0 = 0.12, b0 = 6.;
  printf("\n== 7. Coulomb + Woods-Saxon, alpha on fixed Au (v = %.2f c, b = %.0f fm, grazing) ==\n", v0, b0);
  event* ev = NewEvent(kTRUE);
  ev->AddParticle(new EMparticle(4., 2., TVector3(-1000., b0, 0.), TVector3(v0, 0., 0.)));
  ev->AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), kTRUE));
  const Double_t kappa = kCoulombK*2.*79.;
  // closest approach
  Double_t rmin = 1e300;
  while(kTRUE){
    TVector3 r = ev->getParticle(0)->GetX();
    rmin = TMath::Min(rmin, r.Mag());
    if(r.Dot(ev->getParticle(0)->GetV()) > 0. && r.Mag() > 11.) break;
    ev->DeriveDT(0.05);
  }
  Info("closest approach (fm)", rmin);
  const Int_t nd = 8;
  const Double_t Ds[nd] = {11., 12., 14., 17., 20., 30., 100., 1000.};
  printf("  %6s %10s | %10s %10s | %10s %10s | %10s %10s | %10s\n", "D(fm)", "U_WS/U_C", "Coul.only", "Born",
         "pair", "est", "auto", "est", "brute est");
  Double_t worst = 0., eCoul[nd];
  for(Int_t k=0; k<nd; k++){
    IntegrateTo(ev, Ds[k]);
    particle* p1 = ev->getParticle(0);
    particle* p2 = ev->getParticle(1);
    asymptoteN bru(ev, asymptoteN::kBrute);
    std::vector<TVector3> vref = VInf(bru);
    std::vector<Double_t> eref = Errs(bru);
    asymptoteN born(ev, asymptoteN::kBorn), pair(ev, asymptoteN::kPairwise), aut(ev, asymptoteN::kAuto);
    asymptote coul(p1, p2, kappa);   // ignores the nuclear force
    eCoul[k] = (coul.GetV1Out() - vref[0]).Mag()/v0;
    const Double_t ratio = TMath::Max(RatioVs(pair, vref, eref), RatioVs(aut, vref, eref));
    worst = TMath::Max(worst, TMath::Max(ratio, RatioVs(born, vref, eref)));
    printf("  %6.1f %10.2e | %10.2e %10.2e | %10.2e %10.2e | %10.2e %10.2e | %10.2e\n", p1->GetX().Mag(),
           gWS.Potential(p1, p2)/gCoulomb.Potential(p1, p2), eCoul[k], ErrVs(born, vref)/v0,
           ErrVs(pair, vref)/v0, pair.GetMaxError()/v0, ErrVs(aut, vref)/v0, aut.GetMaxError()/v0, bru.GetMaxError()/v0);
    if(k==0) Info("pairwise run time with the quadrature (ms)", 1e3*pair.GetRunTime());
  }
  printf("  (errors / v0 against the brute force ; Coul.only : exact two-body asymptote ignoring the nuclear force)\n");
  Check("Coulomb-only asymptote error at 11 fm (nuclear force matters) > 1e-6", 1e-6/eCoul[0], 1.);
  Check("Coulomb-only asymptote error at 30 fm (outside its range)", eCoul[5], 1e-11);
  Check("pairwise / auto / Born : max actual / estimated error", worst, 1.);
  gWorstRatio = TMath::Max(gWorstRatio, worst);
  delete ev;
}

////////////////////////////////////////////
// 8. relativistic mode
////////////////////////////////////////////
void Relativistic(){
  const Double_t v0 = 0.3;
  printf("\n== 8. Relativistic mode, alpha on fixed Au (v = %.1f c, b = 20 fm) : first order vs relativistic brute force ==\n", v0);
  event* ev = NewEvent();
  ev->SetRelativistic(kTRUE);
  ev->AddParticle(new EMparticle(4., 2., TVector3(-1000., 20., 0.), TVector3(v0, 0., 0.)));
  ev->AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), kTRUE));
  const Int_t nd = 5;
  const Double_t Ds[nd] = {200., 500., 1000., 3000., 10000.};
  Double_t eB[nd], worst = 0.;
  printf("  %7s %10s %10s %10s %7s %12s\n", "D(fm)", "raw", "Born(rel)", "est", "ratio", "nonrel exact");
  for(Int_t k=0; k<nd; k++){
    IntegrateTo(ev, Ds[k]);
    asymptoteN bru(ev, asymptoteN::kBrute);
    std::vector<TVector3> vref = VInf(bru);
    std::vector<Double_t> eref = Errs(bru);
    asymptoteN a(ev, asymptoteN::kPairwise);   // falls back to first order, flag RELATIVISTIC
    asymptote nr(ev->getParticle(0), ev->getParticle(1), kCoulombK*2.*79.);
    eB[k] = ErrVs(a, vref)/v0;
    const Double_t r = RatioVs(a, vref, eref);
    worst = TMath::Max(worst, r);
    printf("  %7.0f %10.2e %10.2e %10.2e %7.3f %12.2e\n", Ds[k], (ev->getParticle(0)->GetV() - vref[0]).Mag()/v0,
           eB[k], a.GetMaxError()/v0, r, (nr.GetV1Out() - vref[0]).Mag()/v0);
    if(k==0) CheckTrue("flag RELATIVISTIC set", a.GetFlags() & asymptoteN::kRelativistic);
  }
  CheckRange("order of the relativistic first-order error (last interval)", Order(eB[nd-2], eB[nd-1], Ds[nd-2], Ds[nd-1]), 1.8, 2.2);
  Check("relativistic first order : max actual / estimated error", worst, 1.);
  delete ev;
}

////////////////////////////////////////////
// 9. special states and guarantees
////////////////////////////////////////////
void Special(){
  printf("\n== 9. Special states ==\n");
  const Double_t kappa = kCoulombK*2.*79.;
  // two bodies from the APPROACHING start state, b = 50 and head-on
  const Double_t bs[2] = {50., 0.};
  for(Int_t k=0; k<2; k++){
    event* ev = NewEvent();
    ev->AddParticle(new EMparticle(4., 2., TVector3(-1000., bs[k], 0.), TVector3(0.05, 0., 0.)));
    ev->AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), kTRUE));
    asymptote ex(ev->getParticle(0), ev->getParticle(1), kappa);
    asymptoteN p(ev, asymptoteN::kPairwise);
    asymptoteN g; g.SetAnalytic(kFALSE); g.Compute(ev, asymptoteN::kPairwise);
    asymptoteN pin(ev, asymptoteN::kPairwise, asymptoteN::kIncoming);
    TString t = Form("b = %.0f, approaching start : ", bs[k]);
    Check((t + "pairwise closed form vs exact (/ v0)").Data(), (p.GetVInf(0) - ex.GetV1Out()).Mag()/0.05, 1e-14);
    Check((t + "pairwise orbit quadrature vs exact (/ v0)").Data(), (g.GetVInf(0) - ex.GetV1Out()).Mag()/0.05, 1e-12);
    Check((t + "incoming pairwise vs exact (/ v0)").Data(), (pin.GetVInf(0) - ex.GetV1In()).Mag()/0.05, 1e-14);
    CheckTrue((t + "flag APPROACHING, not separating").Data(), (p.GetFlags() & asymptoteN::kApproaching) && !p.IsSeparating());
    if(k==1){
      asymptoteN b(ev, asymptoteN::kBorn);
      CheckTrue((t + "Born flagged SINGULAR (straight line through the centre)").Data(), b.GetFlags() & asymptoteN::kSingular);
    }
    delete ev;
  }
  // the caller's event is not modified (state, step count, time)
  {
    std::vector<spec> s = CaseFree4();
    event* ev = Build(s);
    ev->DeriveDTN(10, 5.);
    std::vector<TVector3> x0, v0;
    for(Int_t i=0; i<ev->getNParticle(); i++){ x0.push_back(ev->getParticle(i)->GetX()); v0.push_back(ev->getParticle(i)->GetV()); }
    const Long_t n0 = ev->GetDeriveN();
    const Double_t t0 = ev->GetTime();
    asymptoteN a1(ev, asymptoteN::kBorn), a2(ev, asymptoteN::kPairwise), a3(ev, asymptoteN::kAuto),
               a4(ev, asymptoteN::kBrute), a5(ev, asymptoteN::kBrute, asymptoteN::kIncoming);
    Bool_t same = (ev->GetDeriveN()==n0) && (ev->GetTime()==t0);
    for(Int_t i=0; i<ev->getNParticle(); i++) same &= (ev->getParticle(i)->GetX()==x0[i]) && (ev->getParticle(i)->GetV()==v0[i]);
    CheckTrue("the caller's event is unchanged by every method (x, v, steps, time)", same);
    delete ev;
  }
  // no force : v_inf = v ; no moving particle / no event : invalid
  {
    event* ev = new event();
    ev->AddParticle(new EMparticle(4., 2., TVector3(-10., 1., 0.), TVector3(0.05, 0., 0.)));
    ev->AddParticle(new EMparticle(4., 2., TVector3(10., 0., 0.), TVector3(-0.05, 0., 0.)));
    asymptoteN a(ev, asymptoteN::kPairwise), b(ev, asymptoteN::kBrute);
    CheckTrue("no force : v_inf = v exactly (pairwise and brute force)",
              a.GetVInf(0)==ev->getParticle(0)->GetV() && b.GetVInf(1)==ev->getParticle(1)->GetV() && a.IsReliable());
    delete ev;
    event* e2 = new event();
    e2->AddForce(&gCoulomb);
    e2->AddParticle(new EMparticle(4., 2., TVector3(), TVector3(), kTRUE));
    asymptoteN c(e2);
    asymptoteN d(nullptr);
    CheckTrue("only fixed particles, or no event : INVALID, Compute returns kFALSE",
              !c.IsValid() && !d.IsValid() && !asymptoteN().Compute(nullptr));
    delete e2;
  }
  // radial, attractive, receding with E > 0 : asymptote gives none, the pair has one
  {
    event* ev = NewEvent();
    ev->AddParticle(new EMparticle(1.,  1., TVector3(50., 0., 0.),  TVector3(0.02, 0., 0.)));
    ev->AddParticle(new EMparticle(1., -1., TVector3(-50., 0., 0.), TVector3(-0.02, 0., 0.)));
    asymptoteN a(ev, asymptoteN::kPairwise), b(ev, asymptoteN::kBrute);
    const Double_t E = 0.5*0.5*0.04*0.04 - kCoulombK/100.;
    const Double_t vexp = 0.5*TMath::Sqrt(2.*E/0.5);
    Check("radial attractive receding pair : |v_inf| vs sqrt(2E/mu)/2 (relative)", TMath::Abs(a.GetVInf(0).X() - vexp)/vexp, 1e-13);
    Check("radial attractive receding pair : pairwise vs brute (/ |v_inf|)", ErrVs(a, VInf(b))/vexp, 1e-10);
    CheckTrue("radial attractive receding pair : reliable", a.IsReliable());
    delete ev;
  }
  // a slow, nearly parabolic pair is flagged
  {
    event* ev = NewEvent();
    const Double_t vpar = TMath::Sqrt(2.*kCoulombK*1.*1./(0.5*100.)); // relative speed for E = 0 at 100 fm (+1, -1)
    ev->AddParticle(new EMparticle(1.,  1., TVector3(50., 0., 0.),  TVector3(0.5*vpar*(1.+1e-14), 0., 0.)));
    ev->AddParticle(new EMparticle(1., -1., TVector3(-50., 0., 0.), TVector3(-0.5*vpar*(1.+1e-14), 0., 0.)));
    asymptoteN a(ev, asymptoteN::kPairwise);
    printf("  nearly parabolic pair : flags %s\n", a.GetFlagNames(a.GetFlags()).Data());
    CheckTrue("nearly parabolic pair (E ~ 1e-14 KE) : MARGINAL, not reliable", (a.GetFlags() & asymptoteN::kMarginal) && !a.IsReliable());
    delete ev;
  }
  // termination criterion in a loop, as a user would write it
  {
    std::vector<spec> s = CaseFixed3();
    event* ev = Build(s);
    Long_t n = 0;
    Bool_t sep = kFALSE;
    while(!(sep && asymptoteN::RemainingDeflection(ev, &sep) < 1e-4)){
      for(Int_t k=0; k<20; k++) ev->DeriveDT(5.);
      n += 20;
      asymptoteN::RemainingDeflection(ev, &sep);
    }
    asymptoteN ref(ev, asymptoteN::kBrute);
    Double_t act = 0.;
    for(Int_t i=1; i<3; i++) act = TMath::Max(act, asymptote::Angle(ev->getParticle(i)->GetV(), ref.GetVInf(i)));
    printf("  stop when RemainingDeflection < 1e-4 rad : after %ld steps (min distance %.0f fm), true remaining %.4e rad\n",
           n, MinPairDistance(ev), act);
    CheckRange("true remaining deflection at that stop / 1e-4", act/1e-4, 0., 1.05);
    asymptoteN b(ev, asymptoteN::kBorn);
    Info("cost of one RemainingDeflection there (microseconds, Coulomb)", 1e6*b.GetRunTime());
    delete ev;
  }
}

int main(int argc, char** argv){
  TwoBody("Two bodies, fixed target", kTRUE, 197., 79., 50.);
  TwoBody("Two bodies, free target", kFALSE, 12., 6., 20.);
  NBody("Three bodies : fixed Au + two alphas", CaseFixed3());
  NBody("Three bodies, all free : C12 + two alphas", CaseFree3());
  NBody("Four bodies, all free : Au + two alphas + proton", CaseFree4());
  Residual("fixed Au + two alphas", CaseFixed3());
  Residual("four free bodies", CaseFree4());
  Incoming("three free bodies (C12 + two alphas)", CaseFree3());
  Incoming("fixed Au + two alphas", CaseFixed3());
  Bound();
  Nuclear2();
  {
    // three bodies with the nuclear force : fixed Au + two grazing alphas at 0.12 c
    std::vector<spec> s;
    s.push_back({197., 79., TVector3(), TVector3(), kTRUE});
    s.push_back({4., 2., TVector3(-1000., 6., 0.), TVector3(0.12, 0., 0.), kFALSE});
    s.push_back({4., 2., TVector3(4., -1100., 3.), TVector3(0., 0.12, 0.), kFALSE});
    NBody("Three bodies with Coulomb + Woods-Saxon (fixed Au + two grazing alphas at 0.12 c)", s, kTRUE, 15., 0.12);
  }
  Relativistic();
  Special();

  printf("\n  largest actual / estimated error :\n");
  printf("    two bodies, all methods (sections 1, 7)       : %.3f  (Born : estimate = |Born - exact pair| by construction)\n", gWorstRatio);
  printf("    3-4 bodies, pairwise and auto (section 2)     : %.3f  (safety factor %.1f)\n", gWorstPairN, 1./gWorstPairN);
  printf("    3-4 bodies, Born (section 2)                  : %.3f  (safety factor %.2f)\n", gWorstBornN, 1./gWorstBornN);
  std::cout<<std::endl<<"== "<<(gCheck-gFail)<<" / "<<gCheck<<" checks passed =="<<std::endl;
  return gFail==0 ? 0 : 1;
}
