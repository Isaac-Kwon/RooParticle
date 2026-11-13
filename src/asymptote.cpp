#include "asymptote.hpp"

#include "iostream"

#include "TMath.h"

#include "particle.hpp"

// Two states closer than this (relative to the scale) are treated as equal.
static const Double_t kHeadOnTol   = 1e-12; // |r x v| <= tol |r||v|      : L = 0
static const Double_t kParabolaTol = 1e-13; // |E| <= tol (KE + |PE|)     : E = 0
static const Double_t kAtRestTol   = 1e-9;  // |v_self| <= tol |v_other|  : at rest

asymptote::asymptote(particle* p1, particle* p2, Double_t kappa_){
  Compute(p1, p2, kappa_);
}

asymptote::asymptote(TVector3 x1, TVector3 v1, Double_t m1_, Bool_t fixed1_,
                     TVector3 x2, TVector3 v2, Double_t m2_, Bool_t fixed2_,
                     Double_t kappa_){
  Compute(x1, v1, m1_, fixed1_, x2, v2, m2_, fixed2_, kappa_);
}

asymptote::asymptote(TVector3 r_, TVector3 v_, Double_t mu_, Double_t kappa_){
  Compute(r_, v_, mu_, kappa_);
}

void asymptote::Reset(){
  st = kInvalid;
  fixedCentre = fixed1 = fixed2 = kFALSE;
  m1 = m2 = mu = kappa = 0.;
  r = v = vcm = TVector3();
  energy = 0.;
  lvec = avec = evec = phat = TVector3();
  ecc = rmin = 0.;
  vinf = 0.;
  vinRel = voutRel = TVector3();
  thetaCM = bimp = 0.;
  v1in = v1out = v2in = v2out = TVector3();
  thetaLab1 = thetaLab2 = 0.;
}

Bool_t asymptote::Compute(particle* p1, particle* p2, Double_t kappa_){
  if(p1==nullptr || p2==nullptr || p1==p2){
    Reset();
    return kFALSE;
  }
  return Compute(p1->GetX(), p1->GetV(), p1->GetM0(), p1->IsInvincible(),
                 p2->GetX(), p2->GetV(), p2->GetM0(), p2->IsInvincible(),
                 kappa_);
}

Bool_t asymptote::Compute(TVector3 x1, TVector3 v1, Double_t m1_, Bool_t fixed1_,
                          TVector3 x2, TVector3 v2, Double_t m2_, Bool_t fixed2_,
                          Double_t kappa_){
  Reset();
  fixed1 = fixed1_;
  fixed2 = fixed2_;
  fixedCentre = (fixed1 != fixed2);
  m1 = m1_;
  m2 = m2_;
  kappa = kappa_;
  if(fixed1 && fixed2) return kFALSE;
  if((!fixed1 && !(m1>0.)) || (!fixed2 && !(m2>0.))) return kFALSE;

  // a fixed centre is at rest whatever its stored velocity is
  if(fixed1) v1 = TVector3();
  if(fixed2) v2 = TVector3();

  r = x1 - x2;
  v = v1 - v2;
  if(fixed2)      { mu = m1; vcm = TVector3(); }
  else if(fixed1) { mu = m2; vcm = TVector3(); }
  else            { mu = m1*m2/(m1+m2); vcm = (v1*m1 + v2*m2)*(1./(m1+m2)); }

  Solve();
  SetLab();
  return HasAsymptote();
}

Bool_t asymptote::Compute(TVector3 r_, TVector3 v_, Double_t mu_, Double_t kappa_){
  // a body of mass mu around a fixed particle 2 at the origin
  return Compute(r_, v_, mu_, kFALSE, TVector3(), TVector3(), 0., kTRUE, kappa_);
}

void asymptote::Solve(){
  st = kInvalid;
  const Double_t rmag = r.Mag();
  if(!(rmag>0.) || !(mu>0.)) return;

  const Double_t v2   = v.Mag2();
  const Double_t vmag = TMath::Sqrt(v2);
  const TVector3 rhat = r*(1./rmag);
  const TVector3 h    = r.Cross(v);           // specific angular momentum
  const Double_t hmag = h.Mag();
  const Double_t ke   = 0.5*mu*v2;
  const Double_t pe   = kappa/rmag;

  energy = ke + pe;
  lvec   = h*mu;
  avec   = v.Cross(h)*(mu*mu) + rhat*(mu*kappa);

  // ---- no interaction : straight line
  if(kappa==0.){
    if(!(vmag>0.)){ st = kParabolic; rmin = rmag; return; } // at rest forever
    st      = kFree;
    vinf    = vmag;
    vinRel  = v;
    voutRel = v;
    thetaCM = 0.;
    bimp    = hmag/vmag;
    rmin    = bimp;
    TVector3 perp = r - v*(r.Dot(v)/v2);    // centre -> closest point of the line
    phat = (perp.Mag()>0.) ? perp.Unit() : TVector3();
    return;
  }

  // gm = -kappa/mu (> 0 attractive). e^2 - 1 = 2 eps h^2 / gm^2, eps = E/mu,
  // computed directly to avoid the cancellation of e^2 - 1 near e = 1.
  const Double_t gm   = -kappa/mu;
  const Double_t agm  = TMath::Abs(gm);
  const Double_t eps  = energy/mu;
  Double_t e2m1 = 2.*eps*hmag*hmag/(gm*gm);
  if(e2m1 < -1.) e2m1 = -1.;                  // round-off of a circular orbit
  ecc  = TMath::Sqrt(1. + e2m1);
  evec = avec*(1./(mu*TMath::Abs(kappa)));
  phat = (avec.Mag2()>0.) ? avec.Unit() : TVector3();

  const Bool_t radial = (hmag <= kHeadOnTol*rmag*vmag);

  // ---- bound and parabolic orbits (attractive only)
  if(TMath::Abs(energy) <= kParabolaTol*(ke + TMath::Abs(pe))){
    st   = radial ? kRadialAttractive : kParabolic;
    rmin = radial ? 0. : hmag*hmag/(2.*agm);
    return;
  }
  if(energy < 0.){
    st   = radial ? kRadialAttractive : kBound;
    rmin = radial ? 0. : (hmag*hmag/agm)/(1.+ecc);
    return;
  }

  // ---- unbound
  vinf = TMath::Sqrt(2.*eps);
  const Double_t a = agm/(vinf*vinf);         // semi-axis of the hyperbola

  if(radial){
    if(kappa < 0.){                            // falls into the centre
      st   = kRadialAttractive;
      vinf = 0.;
      rmin = 0.;
      return;
    }
    // head-on repulsive : comes in and goes out along the line of r
    st      = kHeadOn;
    phat    = rhat;
    voutRel = rhat*vinf;
    vinRel  = rhat*(-vinf);
    thetaCM = TMath::Pi();
    bimp    = 0.;
    rmin    = kappa/energy;                    // = 2a
    return;
  }

  st = kUnbound;
  // Orbital frame : phat (centre -> periapsis), qhat (direction of motion at periapsis).
  // Asymptotic directions, sigma = +1 repulsive, -1 attractive, s = sqrt(e^2-1):
  //   out = ( sigma phat + s qhat) / e,  in = (-sigma phat + s qhat) / e.
  const TVector3 nhat  = h*(1./hmag);
  const TVector3 qhat  = nhat.Cross(phat);
  const Double_t sigma = (kappa > 0.) ? 1. : -1.;
  const Double_t s     = TMath::Sqrt(e2m1);
  voutRel = (phat*( sigma) + qhat*s)*(vinf/ecc);
  vinRel  = (phat*(-sigma) + qhat*s)*(vinf/ecc);
  thetaCM = 2.*TMath::ATan2(1., s);            // sin(theta/2) = 1/e
  bimp    = hmag/vinf;
  if(kappa > 0.) rmin = a*(1.+ecc);
  else           rmin = (hmag*hmag/agm)/(1.+ecc);
}

static Double_t LabAngleOf(const TVector3 & vin, const TVector3 & vout, const TVector3 & vinOther){
  if(vin.Mag() > kAtRestTol*vinOther.Mag() && vin.Mag2()>0.) return asymptote::Angle(vin, vout);
  if(vinOther.Mag2()>0.)                                       return asymptote::Angle(vinOther, vout);
  return 0.;
}

void asymptote::SetLab(){
  if(!HasAsymptote()) return;
  // velocity of each particle = Vcm +- (share of the relative velocity)
  Double_t c1, c2;
  if(fixed2)      { c1 = 1.; c2 = 0.; }
  else if(fixed1) { c1 = 0.; c2 = 1.; }
  else            { c1 = m2/(m1+m2); c2 = m1/(m1+m2); }

  v1in  = vcm + vinRel*c1;
  v1out = vcm + voutRel*c1;
  v2in  = vcm - vinRel*c2;
  v2out = vcm - voutRel*c2;
  if(fixed1){ v1in = v1out = TVector3(); }
  if(fixed2){ v2in = v2out = TVector3(); }

  thetaLab1 = fixed1 ? 0. : LabAngleOf(v1in, v1out, v2in);
  thetaLab2 = fixed2 ? 0. : LabAngleOf(v2in, v2out, v1in);
}

TString asymptote::GetStatusName() const{
  switch(st){
    case kUnbound:          return TString("UNBOUND");
    case kHeadOn:           return TString("HEADON");
    case kFree:             return TString("FREE");
    case kRadialAttractive: return TString("RADIAL_ATTRACTIVE");
    case kParabolic:        return TString("PARABOLIC");
    case kBound:            return TString("BOUND");
    case kInvalid:          return TString("INVALID");
  }
  return TString("");
}

TString asymptote::Print(Bool_t mute) const{
  TString out;
  out += Form("asymptote : %s (%s)\n", GetStatusName().Data(), fixedCentre ? "fixed centre" : "two free bodies");
  out += Form("  kappa %.6e  mu %.6e  E %.9e  |L| %.9e  e %.12f\n", kappa, mu, energy, lvec.Mag(), ecc);
  out += Form("  rmin %.9f  b %.9f  vinf %.12e  thetaCM %.12f deg\n",
              rmin, bimp, vinf, thetaCM*TMath::RadToDeg());
  out += Form("  vout(rel) (%.12e, %.12e, %.12e)\n", voutRel.X(), voutRel.Y(), voutRel.Z());
  out += Form("  v1out (%.12e, %.12e, %.12e)  theta1 %.9f deg\n", v1out.X(), v1out.Y(), v1out.Z(), thetaLab1*TMath::RadToDeg());
  out += Form("  v2out (%.12e, %.12e, %.12e)  theta2 %.9f deg\n", v2out.X(), v2out.Y(), v2out.Z(), thetaLab2*TMath::RadToDeg());
  if(!mute) std::cout<<out;
  return out;
}

// ---- closed-form relations

Double_t asymptote::RutherfordAngle(Double_t b, Double_t vinf, Double_t mu, Double_t kappa){
  return 2.*TMath::ATan2(kappa, mu*vinf*vinf*TMath::Abs(b));
}

Double_t asymptote::RutherfordImpactParameter(Double_t theta, Double_t vinf, Double_t mu, Double_t kappa){
  return TMath::Abs(kappa)/(mu*vinf*vinf*TMath::Tan(0.5*TMath::Abs(theta)));
}

Double_t asymptote::RutherfordClosestApproach(Double_t b, Double_t vinf, Double_t mu, Double_t kappa){
  const Double_t a  = kappa/(mu*vinf*vinf);
  const Double_t sq = TMath::Sqrt(a*a + b*b);
  if(a >= 0.) return a + sq;
  return b*b/(sq - a);                         // = sq + a, without cancellation
}

Double_t asymptote::HeadOnClosestApproach(Double_t vinf, Double_t mu, Double_t kappa){
  return 2.*kappa/(mu*vinf*vinf);
}

Double_t asymptote::LabAngle(Double_t thetaCM, Double_t m1, Double_t m2){
  return TMath::ATan2(TMath::Sin(thetaCM), TMath::Cos(thetaCM) + m1/m2);
}

Double_t asymptote::RecoilAngle(Double_t thetaCM){
  return 0.5*(TMath::Pi() - thetaCM);
}

Double_t asymptote::Angle(const TVector3 & a, const TVector3 & b){
  return TMath::ATan2(a.Cross(b).Mag(), a.Dot(b));
}
