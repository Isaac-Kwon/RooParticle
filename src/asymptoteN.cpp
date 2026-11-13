#include "asymptoteN.hpp"

#include "iostream"
#include "map"
#include "vector"
#include <chrono>
#include <cmath>
#include <typeinfo>

#include "TMath.h"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "asymptote.hpp"

namespace {

typedef asymptoteN::pstate pstate;

// relative round-off floor of every error estimate
const Double_t kRoundOff = 4e-15;
// a non-1/r force below this fraction of the pair potential and force, on a
// receding pair, is taken as finished (short-range forces decay faster than 1/r^2)
const Double_t kShortRange = 1e-17;

////////////////////////////////////////////////////////////////////////////
// quadrature
////////////////////////////////////////////////////////////////////////////
struct glRule { std::vector<Double_t> x, w; };

// Gauss-Legendre nodes and weights on [-1, 1] (Newton on P_n), cached
const glRule & GaussLegendre(Int_t n){
  static std::map<Int_t, glRule> cache;
  std::map<Int_t, glRule>::iterator it = cache.find(n);
  if(it!=cache.end()) return it->second;
  glRule g;
  g.x.resize(n);
  g.w.resize(n);
  for(Int_t i=0; i<n; i++){
    Double_t z = std::cos(TMath::Pi()*(i+0.75)/(n+0.5));
    Double_t pp = 1.;
    for(Int_t iter=0; iter<100; iter++){
      Double_t p1 = 1., p2 = 0.;
      for(Int_t j=1; j<=n; j++){
        Double_t p3 = p2;
        p2 = p1;
        p1 = ((2.*j-1.)*z*p2 - (j-1.)*p3)/j;
      }
      pp = n*(z*p1 - p2)/(z*z - 1.);
      Double_t z1 = z;
      z = z1 - p1/pp;
      if(std::fabs(z-z1) < 1e-15) break;
    }
    // derivative at the final node
    Double_t p1 = 1., p2 = 0.;
    for(Int_t j=1; j<=n; j++){
      Double_t p3 = p2;
      p2 = p1;
      p1 = ((2.*j-1.)*z*p2 - (j-1.)*p3)/j;
    }
    pp = n*(z*p1 - p2)/(z*z - 1.);
    g.x[i] = z;
    g.w[i] = 2./((1.-z*z)*pp*pp);
  }
  return cache[n] = g;
}

// panels on [0, 1]
const Double_t kBrPhi[]  = {0., 0.01, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
const Double_t kBrLow[]  = {0., 1e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 1.};
const Double_t kBrHigh[] = {0., 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.97, 0.99, 0.997, 0.999, 0.9999, 1.};
const Int_t    kNBrPhi   = sizeof(kBrPhi)/sizeof(kBrPhi[0]);
const Int_t    kNBrLow   = sizeof(kBrLow)/sizeof(kBrLow[0]);
const Int_t    kNBrHigh  = sizeof(kBrHigh)/sizeof(kBrHigh[0]);

template<class T, class F>
T Composite(const F & f, const Double_t* br, Int_t nbr, Int_t n){
  const glRule & g = GaussLegendre(n);
  T s = T();
  for(Int_t k=0; k+1<nbr; k++){
    const Double_t a = br[k], b = br[k+1];
    const Double_t hm = 0.5*(b-a), c = 0.5*(a+b);
    for(Int_t i=0; i<n; i++){
      s += f(c + hm*g.x[i]) * (hm*g.w[i]);
    }
  }
  return s;
}

// lower comparison order of the quadrature error estimate
Int_t LowerOrder(Int_t n){ Int_t m = (2*n)/3; return (m<4) ? 4 : m; }

Double_t Norm(Double_t x){ return std::fabs(x); }
Double_t Norm(const TVector3 & v){ return v.Mag(); }

template<class T, class F>
T Panel(const F & f, Double_t a, Double_t b, const glRule & g){
  const Double_t hm = 0.5*(b-a), c = 0.5*(a+b);
  T s = T();
  for(size_t i=0; i<g.x.size(); i++) s += f(c + hm*g.x[i]) * (hm*g.w[i]);
  return s;
}

// Adaptive composite Gauss-Legendre on [0, 1] : each panel is accepted when
// the n- and LowerOrder(n)-point rules agree within tol * (panel width),
// otherwise it is bisected. tol < 0 : relative, |tol| times the size of a
// first non-adaptive pass. The accepted differences are added to err (a
// pessimistic estimate : it is the error of the lower rule).
template<class T, class F>
T Adaptive(const F & f, const Double_t* br, Int_t nbr, Int_t n, Double_t tol, Double_t & err){
  const glRule & gh = GaussLegendre(n);
  const glRule & gl = GaussLegendre(LowerOrder(n));
  if(tol < 0.){
    T s0 = T();
    for(Int_t k=0; k+1<nbr; k++) s0 += Panel<T>(f, br[k], br[k+1], gh);
    tol = -tol*Norm(s0);
    if(!(tol>0.)) tol = 1e-300;
  }
  struct seg { Double_t a, b; Int_t depth; };
  std::vector<seg> stack;
  for(Int_t k=nbr-2; k>=0; k--) stack.push_back({br[k], br[k+1], 0});
  T s = T();
  Int_t npanel = 0;
  while(!stack.empty()){
    seg p = stack.back();
    stack.pop_back();
    const T ih = Panel<T>(f, p.a, p.b, gh);
    const T il = Panel<T>(f, p.a, p.b, gl);
    const Double_t d = Norm(ih - il);
    npanel++;
    // accepted : within tolerance, or at the round-off level of the panel
    // (a 1e-14 kink of the integrand must not exhaust the budget), or too deep
    if(d <= tol*(p.b - p.a) || d <= 1e-12*Norm(ih) || p.depth >= 30 || npanel > 20000){
      s   += ih;
      err += d;
    }else{
      const Double_t m = 0.5*(p.a + p.b);
      stack.push_back({m, p.b, p.depth+1});
      stack.push_back({p.a, m, p.depth+1});
    }
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////
// pair interaction through stand-in particles (mass, charge, position)
////////////////////////////////////////////////////////////////////////////
class pairPot{
public:
  pairPot(const std::vector<force*> & f, Double_t ma_, Double_t qa_, Double_t mb_, Double_t qb_)
    : forces(f), ma(ma_), qa(qa_), mb(mb_), qb(qb_) {}

  // potential at separation r (central force assumed)
  Double_t U(Double_t r) const {
    EMparticle a(ma, qa, TVector3(0., 0., r)), b(mb, qb, TVector3());
    Double_t s = 0.;
    for(size_t k=0; k<forces.size(); k++) s += forces[k]->Potential(&a, &b);
    return s;
  }
  // radial force -dU/dr at separation r
  Double_t Fr(Double_t r) const {
    EMparticle a(ma, qa, TVector3(0., 0., r)), b(mb, qb, TVector3());
    Double_t s = 0.;
    for(size_t k=0; k<forces.size(); k++) s += forces[k]->Evaluate(&a, &b).Z();
    return s;
  }
  // force on a by b, a at d relative to b
  TVector3 F(const TVector3 & d) const {
    EMparticle a(ma, qa, d), b(mb, qb, TVector3());
    TVector3 s;
    for(size_t k=0; k<forces.size(); k++) s += forces[k]->Evaluate(&a, &b);
    return s;
  }
  // force on b by a, a at d relative to b
  TVector3 Fon2(const TVector3 & d) const {
    EMparticle a(ma, qa, d), b(mb, qb, TVector3());
    TVector3 s;
    for(size_t k=0; k<forces.size(); k++) s += forces[k]->Evaluate(&b, &a);
    return s;
  }
  // U(r0) - U(r0 + d) = integral of F_r, without cancellation for small d
  Double_t DeltaU(Double_t r0, Double_t d) const {
    if(std::fabs(d) < 1e-2*r0){
      static const Double_t gx = std::sqrt(0.6);
      const Double_t c = r0 + 0.5*d, hm = 0.5*d;
      return hm*((5./9.)*Fr(c - hm*gx) + (8./9.)*Fr(c) + (5./9.)*Fr(c + hm*gx));
    }
    return U(r0) - U(r0 + d);
  }
private:
  const std::vector<force*> & forces;
  Double_t ma, qa, mb, qb;
};

// Split the forces into the pure 1/r part (unsoftened coulombForce and
// universalForce : U = kappa / r) and the others. Returns kTRUE if there are
// no others.
Bool_t SplitForces(const std::vector<force*> & forces, Double_t ma, Double_t qa,
                   Double_t mb, Double_t qb, Double_t & kappa, std::vector<force*> & others){
  kappa = 0.;
  others.clear();
  for(size_t k=0; k<forces.size(); k++){
    force* f = forces[k];
    if(typeid(*f) == typeid(coulombForce) && static_cast<coulombForce*>(f)->GetSoftening()==0.){
      kappa += static_cast<coulombForce*>(f)->GetK()*qa*qb;
    }else if(typeid(*f) == typeid(universalForce) && static_cast<universalForce*>(f)->GetSoftening()==0.){
      kappa -= static_cast<universalForce*>(f)->GetG()*ma*mb;
    }else{
      others.push_back(f);
    }
  }
  return others.empty();
}
Bool_t PureKappa(const std::vector<force*> & forces, Double_t ma, Double_t qa,
                 Double_t mb, Double_t qb, Double_t & kappa){
  std::vector<force*> others;
  return SplitForces(forces, ma, qa, mb, qb, kappa, others);
}

Double_t PotentialAt(const std::vector<force*> & forces, const pstate & pa, const pstate & pb){
  EMparticle a(pa.m, pa.q, pa.x), b(pb.m, pb.q, pb.x);
  Double_t s = 0.;
  for(size_t k=0; k<forces.size(); k++) s += forces[k]->Potential(&a, &b);
  return s;
}

////////////////////////////////////////////////////////////////////////////
// first order : straight-line impulse on a (r = xa - xb, w = va - vb)
////////////////////////////////////////////////////////////////////////////
TVector3 BornInvR(const TVector3 & r, const TVector3 & w, Double_t kappa, UInt_t & fl){
  if(kappa==0.) return TVector3();
  const Double_t W = w.Mag(), R = r.Mag();
  if(!(W>0.)){ fl |= asymptoteN::kMarginal; return TVector3(); }
  const TVector3 wh = w*(1./W);
  const Double_t s  = r.Dot(wh);
  const TVector3 bv = r - wh*s;
  const Double_t b2 = bv.Mag2();
  Double_t c;
  if(s >= 0.){
    c = 1./(R*(R+s));                 // = (1 - s/R)/b^2, stable
  }else{
    if(b2 <= 1e-20*R*R){ fl |= asymptoteN::kSingular; return TVector3(); }
    c = (1. - s/R)/b2;
  }
  return (bv*c + wh*(1./R))*(kappa/W);
}

TVector3 BornNumeric(const pairPot & pp, const TVector3 & r, const TVector3 & w, Int_t n,
                     Bool_t onB, UInt_t & fl, Double_t & qerr){
  const Double_t W = w.Mag(), R = r.Mag();
  if(!(W>0.)){ fl |= asymptoteN::kMarginal; return TVector3(); }
  Double_t e = 0.;
  const TVector3 wh = w*(1./W);
  const Double_t s  = r.Dot(wh);
  const TVector3 bv = r - wh*s;
  const Double_t b  = bv.Mag();
  auto F = [&](Double_t u)->TVector3{
    TVector3 d = bv + wh*u;
    return onB ? pp.Fon2(d) : pp.F(d);
  };
  TVector3 I;
  if(s >= 0.){
    const Double_t L = R;
    auto f = [&](Double_t y)->TVector3{ Double_t om = 1.-y; return F(s + L*y/om)*(L/(om*om)); };
    I = Adaptive<TVector3>(f, kBrLow, kNBrLow, n, -1e-13, e);
  }else{
    auto f1 = [&](Double_t z)->TVector3{ return F(s*(1.-z))*(-s); };
    I = Adaptive<TVector3>(f1, kBrHigh, kNBrHigh, n, -1e-13, e);
    const Double_t L2 = TMath::Max(b, 1e-3*R);
    auto f2 = [&](Double_t y)->TVector3{ Double_t om = 1.-y; return F(L2*y/om)*(L2/(om*om)); };
    I += Adaptive<TVector3>(f2, kBrLow, kNBrLow, n, -1e-13, e);
  }
  qerr = e/W;
  return I*(1./W);
}

////////////////////////////////////////////////////////////////////////////
// exact isolated two-body problem for a central potential
////////////////////////////////////////////////////////////////////////////
// g(r0 + d) - g(r0),  g(r) = 2 (E - U(r)) / mu - h^2 / r^2  (radial velocity squared)
Double_t GDiff(const pairPot & pp, Double_t mu, Double_t h, Double_t r0, Double_t d){
  const Double_t r = r0 + d;
  return (2./mu)*pp.DeltaU(r0, d) + h*h*d*(2.*r0 + d)/(r0*r0*r*r);
}

// first root inward from R of g(r) = wr2 + GDiff(R, r - R) ; -1 if none
Double_t Turning(const pairPot & pp, Double_t mu, Double_t h, Double_t R, Double_t wr2){
  auto g = [&](Double_t r)->Double_t{ return wr2 + GDiff(pp, mu, h, R, r - R); };
  Double_t hi = R, lo = -1., r = R;
  for(Int_t k=0; k<600; k++){
    r *= 0.95;
    if(g(r) < 0.){ lo = r; break; }
    hi = r;
  }
  if(lo < 0.) return -1.;
  for(Int_t it=0; it<200; it++){
    const Double_t mid = 0.5*(lo + hi);
    if(!(mid>lo && mid<hi)) break;
    if(g(mid) < 0.) lo = mid; else hi = mid;
  }
  return hi;
}

struct pairSol {
  TVector3 wout;
  Double_t rho  = 0.;   // distance for the strength eta (closest approach ahead, or now)
  Double_t qerr = 0.;   // quadrature error of |wout| (c)
  UInt_t   flags = 0;
  Bool_t   ok = kFALSE;
};

pairSol SolveCentral(const pairPot & pp, const TVector3 & r, const TVector3 & w, Double_t mu, Int_t nq){
  pairSol s;
  const Double_t R = r.Mag(), w2 = w.Mag2(), wm = std::sqrt(w2);
  s.rho = R;
  if(!(R>0.) || !(mu>0.)){ s.flags |= asymptoteN::kInvalid; return s; }
  const TVector3 rhat = r*(1./R);
  const Double_t wr   = w.Dot(rhat);
  const TVector3 hv   = r.Cross(w);
  const Double_t h    = hv.Mag();
  const Double_t U0   = pp.U(R);
  const Double_t E    = 0.5*mu*w2 + U0;
  if(!(E > 1e-12*(0.5*mu*w2 + std::fabs(U0)))){ s.flags |= asymptoteN::kMarginal; return s; }
  const Double_t vinf = std::sqrt(2.*E/mu);

  if(h <= 1e-12*R*wm){
    // radial motion
    if(wr >= 0.){ s.wout = rhat*vinf; s.ok = kTRUE; return s; }
    const Double_t rt = Turning(pp, mu, 0., R, wr*wr);
    if(rt > 0.){ s.wout = rhat*vinf; s.rho = rt; }
    else       { s.wout = rhat*(-vinf); s.rho = 0.; s.flags |= asymptoteN::kSingular; }
    s.ok = kTRUE;
    return s;
  }

  const TVector3 qhat = hv.Cross(rhat)*(1./h);  // direction of the transverse motion
  const Bool_t direct = (wr >= 0.3*wm);
  Double_t rmin = R;
  if(!direct){
    rmin = Turning(pp, mu, h, R, wr*wr);
    if(!(rmin>0.)){ s.flags |= asymptoteN::kSingular; return s; }
  }
  // angle swept from r0 to infinity ; r = r0 / (1 - y^2)
  Double_t qe = 0.;
  const Double_t kTolPhi = 1e-14;   // rad
  auto phiInf = [&](Double_t r0, Double_t gbase, Int_t n)->Double_t{
    auto f = [&](Double_t y)->Double_t{
      const Double_t om = 1. - y*y;
      Double_t gg;
      if(om <= 0.) gg = gbase + (2./mu)*pp.U(r0) + h*h/(r0*r0);
      else         gg = gbase + GDiff(pp, mu, h, r0, r0*y*y/om);
      if(!(gg > 0.)) return 0.;
      return 2.*y*(h/r0)/std::sqrt(gg);
    };
    return Adaptive<Double_t>(f, kBrPhi, kNBrPhi, n, kTolPhi, qe);
  };
  // angle swept from rmin to R ; r = rmin + (R - rmin) z^2
  auto phiInner = [&](Int_t n)->Double_t{
    const Double_t D = R - rmin;
    if(!(D>0.)) return 0.;
    auto f = [&](Double_t z)->Double_t{
      const Double_t d = D*z*z, rr = rmin + d;
      const Double_t gg = GDiff(pp, mu, h, rmin, d);
      if(!(gg > 0.)) return 0.;
      return 2.*z*D*h/(rr*rr*std::sqrt(gg));
    };
    return Adaptive<Double_t>(f, kBrLow, kNBrLow, n, kTolPhi, qe);
  };
  auto total = [&](Int_t n)->Double_t{
    if(direct) return phiInf(R, wr*wr, n);
    const Double_t p1 = phiInf(rmin, 0., n);
    if(wr >= 0.)       return p1 - phiInner(n);          // just after the closest approach
    if(wr <= -0.3*wm)  return 2.*p1 - phiInf(R, wr*wr, n); // approaching, far
    return p1 + phiInner(n);                             // approaching, near the turning point
  };
  const Double_t dphi  = total(nq);
  // the estimate is the difference of the two rules (that of the lower one) :
  // pessimistic for the result of the higher one
  s.qerr = vinf*qe;
  s.wout = (rhat*std::cos(dphi) + qhat*std::sin(dphi))*vinf;
  s.rho  = (wr >= 0.) ? R : rmin;
  s.ok   = kTRUE;
  return s;
}

////////////////////////////////////////////////////////////////////////////
// bodies (single particles or bound clusters)
////////////////////////////////////////////////////////////////////////////
struct bodyW {
  std::vector<Int_t> mem;
  Bool_t   fixed  = kFALSE;
  Double_t mdyn   = 0.;   // moving mass
  Double_t mprobe = 0.;   // all rest mass (probe for gravity)
  Double_t q      = 0.;
  TVector3 x, v, p;       // CM position, velocity, momentum (relativistic : sum gamma m v)
  Double_t eint = 0., iso = 0.;
  TVector3 dvPair, dvBorn, dpBorn;
  Double_t cross = 0., qerrPair = 0., qerrBorn = 0., merr = 0., second = 0.;
  UInt_t   flags = 0;
};

void SetBody(bodyW & b, const std::vector<pstate> & ps, Bool_t rel){
  b.fixed = kFALSE;
  b.mdyn = b.mprobe = b.q = 0.;
  b.x = b.v = b.p = TVector3();
  Int_t ifix = -1;
  for(size_t k=0; k<b.mem.size(); k++){
    const pstate & s = ps[b.mem[k]];
    b.mprobe += s.m;
    b.q      += s.q;
    if(s.fixed){ b.fixed = kTRUE; if(ifix<0) ifix = b.mem[k]; continue; }
    b.mdyn += s.m;
    b.x    += s.x*s.m;
    b.v    += s.v*s.m;
    b.p    += rel ? s.v*(s.m/std::sqrt(1. - s.v.Mag2())) : s.v*s.m;
  }
  if(b.fixed){
    b.x = ps[ifix].x;
    b.v = b.p = TVector3();
  }else if(b.mdyn > 0.){
    b.x *= 1./b.mdyn;
    b.v *= 1./b.mdyn;
  }
}

Double_t PairMu(const bodyW & a, const bodyW & b){
  if(a.fixed) return b.mdyn;
  if(b.fixed) return a.mdyn;
  return a.mdyn*b.mdyn/(a.mdyn + b.mdyn);
}

Double_t BodyU(const bodyW & a, const bodyW & b, const std::vector<std::vector<Double_t> > & U,
               const std::vector<pstate> & ps){
  Double_t s = 0.;
  for(size_t i=0; i<a.mem.size(); i++){
    for(size_t j=0; j<b.mem.size(); j++){
      if(ps[a.mem[i]].fixed && ps[b.mem[j]].fixed) continue;
      s += U[a.mem[i]][b.mem[j]];
    }
  }
  return s;
}

void BuildBodies(const std::vector<pstate> & ps, const std::vector<std::vector<Double_t> > & U,
                 Bool_t rel, std::vector<bodyW> & bodies){
  const Int_t n = ps.size();
  bodies.clear();
  for(Int_t i=0; i<n; i++){
    bodyW b;
    b.mem.push_back(i);
    SetBody(b, ps, rel);
    bodies.push_back(b);
  }
  // merge the most bound pair of bodies while one has a negative pair energy
  while(kTRUE){
    Int_t ba = -1, bb = -1;
    Double_t emin = 0.;
    for(size_t a=0; a<bodies.size(); a++){
      for(size_t b=a+1; b<bodies.size(); b++){
        if(bodies[a].fixed && bodies[b].fixed) continue;
        const Double_t u  = BodyU(bodies[a], bodies[b], U, ps);
        const Double_t mu = PairMu(bodies[a], bodies[b]);
        const Double_t ke = 0.5*mu*(bodies[a].v - bodies[b].v).Mag2();
        const Double_t e  = ke + u;
        if(e < -1e-12*(ke + std::fabs(u)) && e < emin){ emin = e; ba = a; bb = b; }
      }
    }
    if(ba < 0) break;
    bodies[ba].mem.insert(bodies[ba].mem.end(), bodies[bb].mem.begin(), bodies[bb].mem.end());
    SetBody(bodies[ba], ps, rel);
    bodies.erase(bodies.begin() + bb);
  }
}

struct pairW {
  Int_t    a = 0, b = 0;
  TVector3 r, w;
  Double_t mu = 0., E = 0., eta = 0.;
  Bool_t   interacting = kFALSE, approaching = kFALSE, havePair = kFALSE, pure = kFALSE;
  Double_t kappa = 0.;
  TVector3 dwPair;          // change of the relative velocity still to come (pairwise or first order)
  TVector3 Ia, Ib;          // first-order impulses on a and b
  Double_t qerrPair = 0., qerrBorn = 0.;
  UInt_t   flags = 0;
};

void ReadState(event* ev, std::vector<pstate> & ps, Bool_t flip){
  const Int_t n = ev->getNParticle();
  ps.resize(n);
  for(Int_t i=0; i<n; i++){
    particle* p = ev->getParticle(i);
    ps[i].m     = p->GetM0();
    ps[i].q     = p->GetQ();
    ps[i].x     = p->GetX();
    ps[i].fixed = p->IsInvincible();
    ps[i].v     = ps[i].fixed ? TVector3() : (flip ? -p->GetV() : p->GetV());
  }
}

// light copy for the integrating methods : stand-in particles (no path
// record), same forces, RK45, same pair symmetry and relativistic dynamics
event* LightCopy(event* ev, Bool_t flip, Double_t rel, Double_t abs){
  event* c = new event();
  for(Int_t k=0; k<ev->getNForce(); k++) c->AddForce(ev->getForce(k));
  for(Int_t i=0; i<ev->getNParticle(); i++){
    particle* p = ev->getParticle(i);
    TVector3 v = p->GetV();
    if(flip && !p->IsInvincible()) v = -v;
    c->AddParticle(new EMparticle(p->GetM0(), p->GetQ(), p->GetX(), v, p->IsInvincible(), kFALSE));
  }
  c->SetPairSymmetric(ev->IsPairSymmetric());
  c->SetRelativistic(ev->IsRelativistic() && ev->GetIntegrator()!=event::kTaylor);
  c->SetIntegrator(event::kRK45);
  c->SetTolerance(rel, abs);
  return c;
}

// time for the closest interacting pair to double its distance
Double_t TimeScale(const std::vector<pstate> & ps, const std::vector<std::vector<Double_t> > & U){
  Double_t t = 0.;
  for(size_t i=0; i<ps.size(); i++){
    for(size_t j=i+1; j<ps.size(); j++){
      if(ps[i].fixed && ps[j].fixed) continue;
      if(U[i][j]==0.) continue;
      const Double_t R = (ps[i].x - ps[j].x).Mag();
      const Double_t W = (ps[i].v - ps[j].v).Mag();
      if(!(W>0.)) continue;
      const Double_t tij = R/W;
      if(t==0. || tij < t) t = tij;
    }
  }
  return t;
}

} // namespace

////////////////////////////////////////////////////////////////////////////
// asymptoteN
////////////////////////////////////////////////////////////////////////////
void asymptoteN::Reset(Int_t n){
  vNow.assign(n, TVector3());
  vBody.assign(n, TVector3());
  vInf.assign(n, TVector3());
  err.assign(n, 0.);
  mass.assign(n, 0.);
  flags.assign(n, 0);
  cluster.assign(n, -1);
  fixed.assign(n, kFALSE);
  clMembers.clear();
  clEnergy.clear();
  clIso.clear();
  clVInf.clear();
  allFlags = 0;
  separating = kTRUE;
  maxDefl = maxErr = vScale = maxEta = etaApp = 0.;
  eNow = eInf = 0.;
  runTime = tInt = 0.;
  nBlocks = 0;
  nSub = 0;
}

void asymptoteN::Analyze(event* ev, const std::vector<pstate> & ps, method m, Bool_t clusterOnly){
  const Int_t n = ps.size();
  std::vector<force*> forces;
  for(Int_t k=0; k<ev->getNForce(); k++) forces.push_back(ev->getForce(k));
  const Bool_t sym = ev->IsPairSymmetric();
  const Bool_t rel = ev->IsRelativistic() && ev->GetIntegrator()!=event::kTaylor;

  // reset the per-analysis outputs
  for(Int_t i=0; i<n; i++){ vInf[i] = TVector3(); vBody[i] = TVector3(); err[i] = 0.; flags[i] = 0; cluster[i] = -1; }
  clMembers.clear(); clEnergy.clear(); clIso.clear(); clVInf.clear();
  allFlags = 0; separating = kTRUE; maxEta = etaApp = 0.;

  Int_t nmoving = 0;
  for(Int_t i=0; i<n; i++){
    mass[i]  = ps[i].m;
    fixed[i] = ps[i].fixed;
    if(!ps[i].fixed){
      nmoving++;
      if(!(ps[i].m > 0.)) allFlags |= kInvalid;
    }
  }
  if(nmoving==0) allFlags |= kInvalid;
  if(allFlags & kInvalid) return;

  // pair potentials at the present positions
  std::vector<std::vector<Double_t> > U(n, std::vector<Double_t>(n, 0.));
  for(Int_t i=0; i<n; i++){
    for(Int_t j=i+1; j<n; j++){
      if(ps[i].fixed && ps[j].fixed) continue;
      U[i][j] = U[j][i] = PotentialAt(forces, ps[i], ps[j]);
    }
  }
  eNow = 0.;
  for(Int_t i=0; i<n; i++){
    if(!ps[i].fixed) eNow += 0.5*ps[i].m*ps[i].v.Mag2();
    for(Int_t j=i+1; j<n; j++) eNow += U[i][j];
  }

  // bodies and clusters
  std::vector<bodyW> bodies;
  BuildBodies(ps, U, rel, bodies);
  const Int_t nb = bodies.size();
  for(Int_t a=0; a<nb; a++){
    bodyW & B = bodies[a];
    if(B.mem.size() < 2) continue;
    // internal energy in the CM frame
    Double_t e = 0.;
    for(size_t i=0; i<B.mem.size(); i++){
      const pstate & s = ps[B.mem[i]];
      if(!s.fixed) e += 0.5*s.m*(s.v - B.v).Mag2();
      for(size_t j=i+1; j<B.mem.size(); j++) e += U[B.mem[i]][B.mem[j]];
    }
    B.eint = e;
    // isolation : tidal / multipole coupling to the other bodies
    Double_t tid = 0.;
    for(Int_t c=0; c<nb; c++){
      if(c==a) continue;
      const bodyW & C = bodies[c];
      if(B.fixed && C.fixed) continue;
      pairPot pp(forces, B.mprobe, B.q, C.mprobe, C.q);
      const Double_t umono = pp.U((B.x - C.x).Mag());
      tid += std::fabs(BodyU(B, C, U, ps) - umono);
    }
    B.iso = (e!=0.) ? tid/std::fabs(e) : 1e300;
    B.flags |= kInCluster;
    if(B.fixed)        B.flags |= kBoundToFixed;
    if(B.iso > isoThr) B.flags |= kNotIsolated;
  }

  if(clusterOnly){
    for(Int_t a=0; a<nb; a++){
      const bodyW & B = bodies[a];
      Int_t cid = -1;
      if(B.mem.size() > 1){
        cid = clMembers.size();
        clMembers.push_back(B.mem);
        clEnergy.push_back(B.eint);
        clIso.push_back(B.iso);
        clVInf.push_back(B.v);
      }
      for(size_t k=0; k<B.mem.size(); k++){
        const Int_t i = B.mem[k];
        vInf[i]    = ps[i].fixed ? TVector3() : B.v;
        vBody[i]   = vInf[i];
        flags[i]   = B.flags;
        cluster[i] = cid;
      }
      allFlags |= B.flags;
    }
    return;
  }

  const Bool_t wantPair = (m!=kBorn) || estimate;

  // pairs of bodies
  std::vector<pairW> pairs;
  for(Int_t a=0; a<nb; a++){
    for(Int_t b=a+1; b<nb; b++){
      bodyW & A = bodies[a];
      bodyW & B = bodies[b];
      if(A.fixed && B.fixed) continue;
      pairW P;
      P.a = a; P.b = b;
      P.r = A.x - B.x;
      P.w = A.v - B.v;
      P.mu = PairMu(A, B);
      const Double_t R = P.r.Mag();
      if(!(R>0.) || !(P.mu>0.)){ allFlags |= kInvalid; continue; }
      pairPot pp(forces, A.mprobe, A.q, B.mprobe, B.q);
      const Double_t U0 = pp.U(R);
      const TVector3 F0 = pp.F(P.r);
      // a neutral cluster has no monopole, but its members still interact
      P.interacting = (U0!=0. || F0.Mag2()!=0.
                       || ((A.mem.size()>1 || B.mem.size()>1) && BodyU(A, B, U, ps)!=0.));
      if(!P.interacting){ pairs.push_back(P); continue; }
      const Bool_t central = F0.Cross(P.r).Mag() <= 1e-10*F0.Mag()*R;
      P.E = 0.5*P.mu*P.w.Mag2() + U0;
      P.approaching = P.r.Dot(P.w) < 0.;
      if(analytic){
        std::vector<force*> others;
        P.pure = SplitForces(forces, A.mprobe, A.q, B.mprobe, B.q, P.kappa, others);
        if(!P.pure && !P.approaching){
          // A receding pair outside the range of its other (short-range)
          // forces : they are finished, the 1/r part is solved in closed form.
          pairPot po(others, A.mprobe, A.q, B.mprobe, B.q);
          const Double_t uo = po.U(R);
          const TVector3 fo = po.F(P.r);
          if(std::fabs(uo) <= kShortRange*std::fabs(U0) && fo.Mag() <= kShortRange*F0.Mag()) P.pure = kTRUE;
        }
      }

      // first order
      if(P.pure){
        P.Ia = BornInvR(P.r, P.w, P.kappa, P.flags);
        P.Ib = -P.Ia;
      }else{
        UInt_t fdummy = 0;
        Double_t qb = 0.;
        P.Ia = BornNumeric(pp, P.r, P.w, nquad, kFALSE, P.flags, P.qerrBorn);
        P.Ib = sym ? -P.Ia : BornNumeric(pp, P.r, P.w, nquad, kTRUE, fdummy, qb);
        P.qerrBorn += qb;
      }

      // exact pair solution
      Double_t rho = R;
      const Bool_t canPair = sym && central;
      if(!canPair) P.flags |= kNoPairwise;
      if(wantPair && canPair){
        Bool_t generic = !P.pure;
        if(P.pure){
          asymptote two(P.r, P.w, P.mu, P.kappa);
          if(two.HasAsymptote()){
            P.dwPair   = two.GetVOutRel() - P.w;
            P.havePair = kTRUE;
            rho = two.IsOutgoing() ? R : two.GetClosestApproach();
          }else if(two.GetStatus()==asymptote::kRadialAttractive && !P.approaching){
            generic = kTRUE;   // radial, attractive, receding : asymptote gives none, the orbit solver does
          }else if(two.GetStatus()==asymptote::kRadialAttractive){
            P.flags |= kSingular;
          }else{
            P.flags |= kMarginal;
          }
        }
        if(generic){
          pairSol s = SolveCentral(pp, P.r, P.w, P.mu, nquad);
          P.flags |= s.flags;
          if(s.ok){
            P.dwPair   = s.wout - P.w;
            P.qerrPair = s.qerr;
            P.havePair = kTRUE;
            rho = s.rho;
          }
        }
      }
      if(!P.havePair){
        // first order in place of the pair solution
        TVector3 dva = A.fixed ? TVector3() : P.Ia*(1./A.mdyn);
        TVector3 dvb = B.fixed ? TVector3() : P.Ib*(1./B.mdyn);
        P.dwPair = dva - dvb;
      }
      // strength
      if(P.E > 1e-12*(0.5*P.mu*P.w.Mag2() + std::fabs(U0))){
        Double_t num;
        if(rho > 0.) num = TMath::Max(std::fabs(pp.U(rho)), rho*std::fabs(pp.Fr(rho)));
        else         num = 1e300;
        P.eta = num/(2.*P.E);
      }else{
        P.eta = 1e300;
        P.flags |= kMarginal;
      }
      if(P.approaching){
        P.flags |= kApproaching;
        if(P.eta > 1e-14) separating = kFALSE;
        etaApp = TMath::Max(etaApp, P.eta);
      }
      maxEta = TMath::Max(maxEta, P.eta);
      pairs.push_back(P);
    }
  }

  // sums over the pairs
  for(size_t k=0; k<pairs.size(); k++){
    const pairW & P = pairs[k];
    if(!P.interacting) continue;
    bodyW & A = bodies[P.a];
    bodyW & B = bodies[P.b];
    if(!A.fixed){
      A.dvPair += P.dwPair*(P.mu/A.mdyn);
      A.dvBorn += P.Ia*(1./A.mdyn);
      A.dpBorn += P.Ia;
      A.flags  |= P.flags;
      A.qerrPair += P.qerrPair*(P.mu/A.mdyn);
      A.qerrBorn += P.qerrBorn/A.mdyn;
    }
    if(!B.fixed){
      B.dvPair -= P.dwPair*(P.mu/B.mdyn);
      B.dvBorn += P.Ib*(1./B.mdyn);
      B.dpBorn += P.Ib;
      B.flags  |= P.flags;
      B.qerrPair += P.qerrPair*(P.mu/B.mdyn);
      B.qerrBorn += P.qerrBorn/B.mdyn;
    }
  }

  // error terms
  for(size_t k=0; k<pairs.size(); k++){
    const pairW & P = pairs[k];
    if(!P.interacting) continue;
    bodyW & A = bodies[P.a];
    bodyW & B = bodies[P.b];
    const TVector3 dVa = A.fixed ? TVector3() : A.dvPair;
    const TVector3 dVb = B.fixed ? TVector3() : B.dvPair;
    // cross term : the other bodies change the relative velocity of this pair
    const Double_t ext  = ((dVa - dVb) - P.dwPair).Mag();
    const Double_t term = (P.eta < 1e299) ? crossC*P.eta*ext : 1e300;
    // second order of a pair without its exact solution
    const Double_t sec  = P.havePair ? 0. : ((P.eta < 1e299) ? P.eta*P.dwPair.Mag() : 1e300);
    const Bool_t strong = (crossC*P.eta > 0.5) && (ext > 1e-12*P.w.Mag());
    if(!A.fixed){ A.cross += term*(P.mu/A.mdyn); A.second += sec*(P.mu/A.mdyn); if(strong) A.flags |= kStrong; }
    if(!B.fixed){ B.cross += term*(P.mu/B.mdyn); B.second += sec*(P.mu/B.mdyn); if(strong) B.flags |= kStrong; }
    // multipole error of clusters : members (moving rigidly) vs point composite
    if(A.mem.size()>1 || B.mem.size()>1){
      TVector3 isum;
      for(size_t i=0; i<A.mem.size(); i++){
        for(size_t j=0; j<B.mem.size(); j++){
          const pstate & si = ps[A.mem[i]];
          const pstate & sj = ps[B.mem[j]];
          if(si.fixed && sj.fixed) continue;
          pairPot pm(forces, si.m, si.q, sj.m, sj.q);
          UInt_t fd = 0;
          Double_t kap = 0.;
          if(analytic && PureKappa(forces, si.m, si.q, sj.m, sj.q, kap)) isum += BornInvR(si.x - sj.x, P.w, kap, fd);
          else{ Double_t qd = 0.; isum += BornNumeric(pm, si.x - sj.x, P.w, nquad, kFALSE, fd, qd); }
        }
      }
      const Double_t dI = (isum - P.Ia).Mag();
      if(!A.fixed) A.merr += dI/A.mdyn;
      if(!B.fixed) B.merr += dI/B.mdyn;
    }
  }

  // results
  const Bool_t useBorn = (m==kBorn) || rel;
  // velocity scale of the event for the round-off floor (relative velocities
  // of a slow body are differences of the larger velocities)
  Double_t vglob = 0.;
  for(Int_t a=0; a<nb; a++){
    if(bodies[a].fixed) continue;
    vglob = TMath::Max(vglob, bodies[a].v.Mag());
    vglob = TMath::Max(vglob, (bodies[a].v + bodies[a].dvPair).Mag());
  }
  for(Int_t a=0; a<nb; a++){
    bodyW & B = bodies[a];
    if(rel) B.flags |= kRelativistic;
    TVector3 vinf;
    Double_t e = 0.;
    if(!B.fixed){
      if(useBorn){
        if(rel){
          const TVector3 p = B.p + B.dpBorn;
          vinf = p*(1./std::sqrt(B.mdyn*B.mdyn + p.Mag2()));
        }else{
          vinf = B.v + B.dvBorn;
        }
        e = B.cross + B.second + B.qerrBorn + B.merr;
        if(wantPair) e += (B.dvBorn - B.dvPair).Mag();
      }else{
        vinf = B.v + B.dvPair;
        e = B.cross + B.second + B.qerrPair + B.merr;
      }
    }
    // round-off floor : a few ulp of the velocities involved
    if(!B.fixed) e += kRoundOff*vglob;
    Int_t cid = -1;
    if(B.mem.size() > 1){
      cid = clMembers.size();
      clMembers.push_back(B.mem);
      clEnergy.push_back(B.eint);
      clIso.push_back(B.iso);
      clVInf.push_back(vinf);
    }
    for(size_t k=0; k<B.mem.size(); k++){
      const Int_t i = B.mem[k];
      vInf[i]    = ps[i].fixed ? TVector3() : vinf;
      vBody[i]   = ps[i].fixed ? TVector3() : B.v;
      err[i]     = ps[i].fixed ? 0. : e;
      flags[i]   = B.flags;
      cluster[i] = cid;
    }
    allFlags |= B.flags;
  }
}

void asymptoteN::RunIntegrating(event* ev, method m){
  const Bool_t flip = (dir==kIncoming);
  event* c = LightCopy(ev, flip, rkRel, rkAbs);
  std::vector<pstate> ps;
  ReadState(c, ps, kFALSE);
  const Int_t n = ps.size();

  std::vector<force*> forces;
  for(Int_t k=0; k<ev->getNForce(); k++) forces.push_back(ev->getForce(k));
  std::vector<std::vector<Double_t> > U(n, std::vector<Double_t>(n, 0.));
  for(Int_t i=0; i<n; i++){
    for(Int_t j=i+1; j<n; j++){
      if(ps[i].fixed && ps[j].fixed) continue;
      U[i][j] = U[j][i] = PotentialAt(forces, ps[i], ps[j]);
    }
  }
  const Double_t T0 = TimeScale(ps, U);

  auto intErr = [&](Double_t scale)->Double_t{
    const Long_t ns = c->GetRK45NSubsteps() + c->GetRK45NRejected();
    return (rkRel*scale + rkAbs)*std::sqrt((Double_t) ns);
  };
  auto scaleOf = [&](const std::vector<TVector3> & v)->Double_t{
    Double_t s = 0.;
    for(Int_t i=0; i<n; i++) if(!fixed[i]) s = TMath::Max(s, v[i].Mag());
    return s;
  };

  if(m==kAuto){
    for(Int_t k=0; ; k++){
      Analyze(c, ps, kPairwise);
      if(allFlags & kInvalid) break;
      const Double_t sc = scaleOf(vInf);
      Double_t emax = 0.;
      for(Int_t i=0; i<n; i++) emax = TMath::Max(emax, err[i]);
      const Bool_t ok = (emax + intErr(sc) <= autoTol*sc)
                        && !(allFlags & (kStrong | kSingular | kMarginal))
                        && !(etaApp > 1e-4);
      if(ok) break;
      if(k>=maxBlocks || !(T0>0.)){ allFlags |= kNotConverged; break; }
      const Double_t dt = T0*TMath::Power(2., k);
      c->DeriveDT(dt);
      tInt += dt;
      nBlocks++;
      ReadState(c, ps, kFALSE);
    }
    const Double_t ie = intErr(scaleOf(vInf));
    for(Int_t i=0; i<n; i++){
      if(fixed[i]) continue;
      err[i] += ie;
      if(allFlags & kNotConverged) flags[i] |= kNotConverged;
    }
  }else{
    // brute force with Richardson extrapolation in 1/t
    std::vector<TVector3> vPrev, rPrev;
    std::vector<Int_t>    sigPrev;
    Double_t tPrev = 0.;
    Bool_t   haveV = kFALSE, haveR = kFALSE, done = kFALSE;
    if(!(T0>0.)){
      Analyze(c, ps, kBrute, kTRUE);
      done = kTRUE;
    }
    for(Int_t k=0; !done; k++){
      const Double_t dt = T0*TMath::Power(2., k);
      c->DeriveDT(dt);
      tInt += dt;
      nBlocks++;
      ReadState(c, ps, kFALSE);
      Analyze(c, ps, kBrute, kTRUE);
      if(allFlags & kInvalid) break;
      const Double_t t = tInt + T0;
      std::vector<TVector3> vk = vInf;
      std::vector<Int_t>    sig = cluster;
      if(haveV && sig!=sigPrev){ haveV = haveR = kFALSE; }   // cluster structure changed : restart
      if(haveV){
        std::vector<TVector3> rk(n);
        for(Int_t i=0; i<n; i++) rk[i] = (vk[i]*t - vPrev[i]*tPrev)*(1./(t - tPrev));
        if(haveR){
          const Double_t sc = scaleOf(rk);
          Double_t dmax = 0.;
          std::vector<Double_t> d(n, 0.);
          for(Int_t i=0; i<n; i++){ if(fixed[i]) continue; d[i] = (rk[i] - rPrev[i]).Mag(); dmax = TMath::Max(dmax, d[i]); }
          const Bool_t conv = (k>=3) && dmax <= bruteTol*sc;
          if(conv || k>=maxBlocks){
            const Double_t ie = intErr(sc);
            for(Int_t i=0; i<n; i++){
              vInf[i] = fixed[i] ? TVector3() : rk[i];
              err[i]  = fixed[i] ? 0. : d[i] + ie;
              if(!conv && !fixed[i]) flags[i] |= kNotConverged;
            }
            for(size_t cc=0; cc<clMembers.size(); cc++) clVInf[cc] = vInf[clMembers[cc][0]];
            if(!conv) allFlags |= kNotConverged;
            done = kTRUE;
          }
        }
        rPrev = rk;
        haveR = kTRUE;
      }
      vPrev = vk;
      sigPrev = sig;
      tPrev = t;
      haveV = kTRUE;
    }
  }
  nSub = c->GetRK45NSubsteps() + c->GetRK45NRejected();
  delete c;
}

void asymptoteN::Finish(event* ev){
  const Int_t n = vInf.size();
  if(dir==kIncoming){
    for(Int_t i=0; i<n; i++){ vInf[i] = -vInf[i]; vBody[i] = -vBody[i]; }
    for(size_t c=0; c<clVInf.size(); c++) clVInf[c] = -clVInf[c];
  }
  for(Int_t i=0; i<n; i++){
    particle* p = ev->getParticle(i);
    vNow[i] = p->IsInvincible() ? TVector3() : p->GetV();
  }
  maxDefl = maxErr = vScale = 0.;
  eInf = 0.;
  for(Int_t i=0; i<n; i++){
    if(fixed[i]) continue;
    maxDefl = TMath::Max(maxDefl, GetDeflection(i));
    maxErr  = TMath::Max(maxErr, err[i]);
    vScale  = TMath::Max(vScale, vInf[i].Mag());
    eInf   += 0.5*mass[i]*vInf[i].Mag2();
  }
  for(size_t c=0; c<clEnergy.size(); c++) eInf += clEnergy[c];
}

Bool_t asymptoteN::Compute(event* ev, method m, direction d){
  auto t0 = std::chrono::steady_clock::now();
  meth = m;
  dir  = d;
  if(ev==nullptr || ev->getNParticle()==0){
    Reset(0);
    allFlags = kInvalid;
    return kFALSE;
  }
  const Int_t n = ev->getNParticle();
  Reset(n);
  std::vector<pstate> ps;
  ReadState(ev, ps, d==kIncoming);
  if(m==kBorn || m==kPairwise){
    Analyze(ev, ps, m);
  }else{
    // velocities of the bodies and the energy of the given state
    Analyze(ev, ps, kPairwise, kTRUE);
    std::vector<TVector3> vb = vBody;
    const Double_t e0 = eNow;
    if(!(allFlags & kInvalid)) RunIntegrating(ev, m);
    vBody = vb;
    eNow  = e0;
  }
  Finish(ev);
  runTime = std::chrono::duration<Double_t>(std::chrono::steady_clock::now() - t0).count();
  return IsValid();
}

Double_t asymptoteN::GetDeflection(Int_t i) const{
  const TVector3 & a = vBody.at(i);
  const TVector3 & b = vInf.at(i);
  if(!(a.Mag2()>0.) || !(b.Mag2()>0.)) return 0.;
  return asymptote::Angle(a, b);
}

Double_t asymptoteN::GetAngularError(Int_t i) const{
  const Double_t v = vInf.at(i).Mag();
  return (v>0.) ? err.at(i)/v : 0.;
}

TVector3 asymptoteN::GetMomentumInf() const{
  TVector3 p;
  for(size_t i=0; i<vInf.size(); i++) if(!fixed[i]) p += vInf[i]*mass[i];
  return p;
}

TString asymptoteN::GetFlagNames(UInt_t f) const{
  static const char* names[] = {"APPROACHING", "STRONG", "CLUSTER", "NOT_ISOLATED", "BOUND_TO_FIXED",
                                "RELATIVISTIC", "NO_PAIRWISE", "MARGINAL", "SINGULAR", "NOT_CONVERGED", "INVALID"};
  TString s;
  for(Int_t k=0; k<11; k++){
    if(f & (1u<<k)){
      if(s.Length()>0) s += "|";
      s += names[k];
    }
  }
  if(s.Length()==0) s = "OK";
  return s;
}

TString asymptoteN::Print(Bool_t mute) const{
  static const char* mnames[] = {"BORN", "PAIRWISE", "AUTO", "BRUTE"};
  TString out;
  out += Form("asymptoteN : %s, %s, flags %s\n", mnames[meth], dir==kOutgoing ? "outgoing" : "incoming",
              GetFlagNames(allFlags).Data());
  out += Form("  max deflection %.3e rad, max error %.3e c (relative %.3e), max eta %.3e, separating %s\n",
              maxDefl, maxErr, GetMaxRelError(), maxEta, separating ? "yes" : "no");
  out += Form("  energy now %.12e, at infinity %.12e ; run time %.3e s, integrated %.3e fm/c, %ld substeps\n",
              eNow, eInf, runTime, tInt, nSub);
  for(size_t i=0; i<vInf.size(); i++){
    out += Form("  %2d %s v_inf (% .12e, % .12e, % .12e) err %.2e defl %.3e cl %d %s\n", (Int_t) i,
                fixed[i] ? "fixed " : "moving", vInf[i].X(), vInf[i].Y(), vInf[i].Z(), err[i],
                GetDeflection(i), cluster[i], GetFlagNames(flags[i]).Data());
  }
  if(!mute) std::cout<<out;
  return out;
}

Double_t asymptoteN::RemainingDeflection(event* ev, Bool_t* sep){
  asymptoteN a;
  a.SetErrorEstimate(kFALSE);
  a.Compute(ev, kBorn);
  if(sep) *sep = a.IsSeparating();
  return a.GetMaxDeflection();
}

Double_t asymptoteN::RemainingChange(event* ev, Bool_t* sep){
  asymptoteN a;
  a.SetErrorEstimate(kFALSE);
  a.Compute(ev, kBorn);
  if(sep) *sep = a.IsSeparating();
  Double_t r = 0.;
  for(Int_t i=0; i<a.GetNParticle(); i++){
    if(a.IsFixed(i)) continue;
    const Double_t v = a.GetVInf(i).Mag();
    if(v>0.) r = TMath::Max(r, a.GetResidual(i).Mag()/v);
  }
  return r;
}
