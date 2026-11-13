#include "iostream"
#include "list"
#include "vector"
#include "stdio.h"
#include <limits>

#include "RPConfig.hpp"

#include "TROOT.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "TVector3.h"

#include "force.hpp"

#include "particle.hpp"
#include "event.hpp"
#include "inspector.hpp"

event::event(){
}

//Copy Constructor
event::event(event &e){
  //Copy all particle into new event
  vector<particle*>::iterator p_;
  for(p_=e.particles.begin(); p_!=e.particles.end(); p_++){
    AddParticle(new particle(**p_));
  }

  //Add all force into new event (not duplicate)
  vector<force*>::iterator ff_;
  for(ff_=e.forces.begin(); ff_!=e.forces.end(); ff_++){
    AddForce(*ff_);
  }

  //copy inspectors : link all particles in each inspectors
  vector<inspector*>::iterator ii_;

  for(ii_=e.inspectors.begin(); ii_!=e.inspectors.end(); ii_++){
    if((*ii_)->IsInspectorP()){
      inspectorP* iis; //source inspector
      iis = (inspectorP*) *ii_;
      Int_t i;
      // found again for every inspector : an index left from the previous
      // inspector would silently link this one to the wrong particles.
      Int_t ip1=-1, ip2=-1;
      for(i=0; i<(Int_t) e.particles.size();i++){
        if((iis->GetParticle1())==e.particles[i]){
          ip1 = i;
        }
        if((iis->GetParticle2())==e.particles[i]){
          ip2 = i;
        }
      }
      if(ip1<0 || ip2<0){
        std::cerr<<"event::event(event&) : inspector ("<<iis->GetMethod()
                 <<") watches a particle that is not in the copied event. It is not copied."<<std::endl;
        continue;
      }
      inspectorP* iin = new inspectorP( // Make new inspector (linking particles)
        getParticle(ip1), getParticle(ip2),
        iis->GetMethod(), iis->GetCriterion());
      iin->SetAbort(iis->IsAbort());
      AddInspector(iin);
    }else if((*ii_)->IsInspectorE()){
      // event inspector : the copy watches the new event, not the source one
      inspectorE* iie = (inspectorE*) *ii_;
      inspectorE* iin = new inspectorE(this, iie->GetMethod(), iie->GetCriterion());
      iin->SetAbort(iie->IsAbort());
      AddInspector(iin);
    }
  }

  //copy other properties
  volume = e.volume;
  nparticle = e.nparticle;
  nforce = e.nforce;
  integ = e.integ;
  imode = e.imode;
  maxInspectSteps = e.maxInspectSteps;
  pairSymmetric = e.pairSymmetric;
  relativistic = e.relativistic;
  rk45Rel = e.rk45Rel;
  rk45Abs = e.rk45Abs;
  // rk45H is not copied : the copy finds its own substep from its first call.
}

event::~event(){
  // The event owns its particles, its inspectors and the pre-derived copy.
  // Forces and the volume are shared with the copies : not deleted here.
  vector<particle*>::iterator p_;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    delete *p_;
  }
  particles.clear();

  vector<inspector*>::iterator i_;
  for(i_=inspectors.begin(); i_!=inspectors.end(); i_++){
    delete *i_;
  }
  inspectors.clear();

  delete pEvent;
  pEvent = nullptr;
}

void event::SetTolerance(Double_t rel, Double_t abs){
  // A zero (or negative) tolerance can never be met : keep the previous one.
  if(!(rel>=0.) || !(abs>=0.) || (rel==0. && abs==0.)){
    std::cerr<<"event::SetTolerance : invalid tolerance (rel "<<rel<<", abs "<<abs
             <<"). Tolerance is not changed."<<std::endl;
    return;
  }
  rk45Rel = rel;
  rk45Abs = abs;
}

void event::CalibrateVelocity(){
  vector<particle*>::iterator p_;
  particle *p;

  for(p_=particles.begin(); p_!=particles.end(); p_++){
    p=*p_;
    if(p->IsInvincible()) continue;
    p->calibrateVelocity(this);
  }
}

Bool_t event::SetIntegrator(TString name){
  name.ToUpper();
  if(name=="TAYLOR"){
    SetIntegrator(kTaylor);
  }else if(name=="SEULER"){
    SetIntegrator(kSymplecticEuler);
  }else if(name=="VERLET"){
    SetIntegrator(kVelocityVerlet);
  }else if(name=="RK4"){
    SetIntegrator(kRK4);
  }else if(name=="RK45"){
    SetIntegrator(kRK45);
  }else{
    std::cerr<<"event::SetIntegrator : No integrator named / "<<name
             <<" (TAYLOR, SEULER, VERLET, RK4, RK45). Integrator is not changed."<<std::endl;
    return kFALSE;
  }
  return kTRUE;
}

TString event::GetIntegratorName(){
  switch(integ){
    case kTaylor:          return TString("TAYLOR");
    case kSymplecticEuler: return TString("SEULER");
    case kVelocityVerlet:  return TString("VERLET");
    case kRK4:             return TString("RK4");
    case kRK45:            return TString("RK45");
  }
  return TString("");
}

TString event::GetStopReasonName(){
  switch(lastStop){
    case kStopNone:      return TString("NONE");
    case kStopInspector: return TString("INSPECTOR");
    case kStopAbort:     return TString("ABORT");
    case kStopMaxSteps:  return TString("MAXSTEPS");
    case kStopInvalid:   return TString("INVALID");
    case kStopFailed:    return TString("FAILED");
  }
  return TString("");
}

void event::DeriveDT(Double_t dt){
  switch(integ){
    case kSymplecticEuler: DeriveSymplecticEuler(dt); break;
    case kVelocityVerlet:  DeriveVelocityVerlet(dt);  break;
    case kRK4:             DeriveRK4(dt);             break;
    case kRK45:            DeriveRK45(dt);            break;
    case kTaylor:
    default:               DeriveTaylor(dt);          break;
  }
}

void event::ComputeForce(std::vector<TVector3> & frc){
  const size_t n = particles.size();
  frc.assign(n, TVector3());

  size_t i, j;
  vector<force*>::iterator ff_;
  if(pairSymmetric){
    // Each unordered pair once : the force on p1 by p2 goes to p1, its opposite to p2.
    // The contributions to one particle are still summed in the order of the
    // particles, so only the sign of the evaluation differs from the other branch.
    for(i=0; i<n; i++){
      particle* p1 = particles[i];
      Bool_t move1 = !p1->IsInvincible();
      for(j=i+1; j<n; j++){
        particle* p2 = particles[j];
        Bool_t move2 = !p2->IsInvincible();
        if(!move1 && !move2) continue;
        for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
          if(move1){
            TVector3 f = (*ff_)->Evaluate(p1, p2);
            frc[i] += f;
            if(move2) frc[j] -= f;
          }else{
            // p1 is fixed : only the force on p2 is needed
            frc[j] += (*ff_)->Evaluate(p2, p1);
          }
        }
      }
    }
  }else{
    // Every ordered pair : for forces without action-reaction symmetry.
    for(i=0; i<n; i++){
      particle* p1 = particles[i];
      if(p1->IsInvincible()) continue;
      for(j=0; j<n; j++){
        if(i==j) continue;
        for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
          frc[i] += (*ff_)->Evaluate(p1, particles[j]);
        }
      }
    }
  }
}

void event::ComputeAcceleration(std::vector<TVector3> & acc){
  ComputeForce(acc);
  size_t i;
  for(i=0; i<acc.size(); i++){
    particle* p = particles[i];
    if(p->IsInvincible()) continue;
    acc[i] = acc[i] * (1./p->GetM(kTRUE));
  }
}

TVector3 event::MomentumOf(particle * p, const TVector3 & v){
  if(!UseRelativistic()) return v;
  return v * (p->GetM0() / TMath::Sqrt(1. - v.Mag2()));
}

TVector3 event::VelocityOf(particle * p, const TVector3 & u){
  if(!UseRelativistic()) return u;
  Double_t m0 = p->GetM0();
  return u * (1. / TMath::Sqrt(m0*m0 + u.Mag2()));
}

void event::GetState(std::vector<TVector3> & x, std::vector<TVector3> & u){
  const size_t n = particles.size();
  x.resize(n);
  u.resize(n);
  size_t i;
  for(i=0; i<n; i++){
    x[i] = particles[i]->GetX();
    u[i] = MomentumOf(particles[i], particles[i]->GetV());
  }
}

void event::SetState(const std::vector<TVector3> & x, const std::vector<TVector3> & u){
  size_t i;
  for(i=0; i<particles.size(); i++){
    particle* p = particles[i];
    if(p->IsInvincible()) continue;
    p->applyMechanic(x[i], VelocityOf(p, u[i]));
  }
}

void event::ComputeDerivative(std::vector<TVector3> & dx, std::vector<TVector3> & du){
  // dx/dt = v ; du/dt = F (relativistic, u = p) or F/m0 (u = v)
  if(UseRelativistic()){
    ComputeForce(du);
  }else{
    ComputeAcceleration(du);
  }
  const size_t n = particles.size();
  dx.resize(n);
  size_t i;
  for(i=0; i<n; i++){
    dx[i] = particles[i]->GetV();
  }
}

void event::FinishStep(Double_t dt){
  vector<particle*>::iterator pp_;
  for(pp_=particles.begin(); pp_!=particles.end(); pp_++){
    (*pp_)->finishStep(dt);
  }
}

void event::DeriveSymplecticEuler(Double_t dt){
  size_t i;
  if(UseRelativistic()){
    // p(t+dt) = p(t) + F(t) dt ; x(t+dt) = x(t) + v(p(t+dt)) dt
    std::vector<TVector3> frc;
    ComputeForce(frc);
    for(i=0; i<particles.size(); i++){
      particle* p = particles[i];
      if(p->IsInvincible()) continue;
      p->applyV(VelocityOf(p, MomentumOf(p, p->GetV()) + frc[i] * dt));
      p->applyDX(p->GetV() * dt);
    }
    FinishStep(dt);
    return;
  }

  // v(t+dt) = v(t) + a(t) dt ; x(t+dt) = x(t) + v(t+dt) dt
  std::vector<TVector3> acc;
  ComputeAcceleration(acc);

  for(i=0; i<particles.size(); i++){
    particle* p = particles[i];
    if(p->IsInvincible()) continue;
    p->applyDV(acc[i] * dt);
    p->applyDX(p->GetV() * dt);
  }
  FinishStep(dt);
}

void event::DeriveVelocityVerlet(Double_t dt){
  // x(t+dt) = x(t) + v(t) dt + a(t) dt^2/2
  // v(t+dt) = v(t) + (a(t) + a(t+dt)) dt/2
  const size_t n = particles.size();
  size_t i;

  // a(t) : reuse the one calculated at the end of the previous step,
  // only if every particle is still where that step left it.
  Bool_t cached = (accCache.size()==n && accCacheX.size()==n);
  for(i=0; cached && i<n; i++){
    if(particles[i]->GetX() != accCacheX[i]) cached = kFALSE;
  }
  if(UseRelativistic()){
    // Kick-drift-kick on the momentum (same scheme as below, as x += v dt + a dt^2/2
    // is x += v(t+dt/2) dt) :
    //   p(t+dt/2) = p(t) + F(t) dt/2 ; x(t+dt) = x(t) + v(p(t+dt/2)) dt
    //   p(t+dt)   = p(t+dt/2) + F(t+dt) dt/2
    // The cache holds forces here ; SetRelativistic() clears it.
    if(!cached){
      ComputeForce(accCache);
    }
    std::vector<TVector3> uh(n);
    for(i=0; i<n; i++){
      particle* p = particles[i];
      if(p->IsInvincible()) continue;
      uh[i] = MomentumOf(p, p->GetV()) + accCache[i] * (0.5*dt);
      TVector3 vh = VelocityOf(p, uh[i]);
      p->applyMechanic(p->GetX() + vh * dt, vh);
    }

    ComputeForce(accCache);

    accCacheX.resize(n);
    for(i=0; i<n; i++){
      particle* p = particles[i];
      accCacheX[i] = p->GetX();
      if(p->IsInvincible()) continue;
      p->applyV(VelocityOf(p, uh[i] + accCache[i] * (0.5*dt)));
    }
    FinishStep(dt);
    return;
  }

  if(!cached){
    ComputeAcceleration(accCache);
  }
  std::vector<TVector3> acc0 = accCache;

  for(i=0; i<n; i++){
    particle* p = particles[i];
    if(p->IsInvincible()) continue;
    p->applyDX(p->GetV() * dt + acc0[i] * (0.5*dt*dt));
  }

  ComputeAcceleration(accCache);

  accCacheX.resize(n);
  for(i=0; i<n; i++){
    particle* p = particles[i];
    accCacheX[i] = p->GetX();
    if(p->IsInvincible()) continue;
    p->applyDV((acc0[i] + accCache[i]) * (0.5*dt));
  }
  FinishStep(dt);
}

void event::DeriveRK4(Double_t dt){
  // Classical 4th order Runge-Kutta on (x, v) of all non-invincible particles.
  // The state is (x, u) : u = v, or u = p in the relativistic mode.
  const size_t n = particles.size();
  size_t i;

  std::vector<TVector3> x0, u0, xs(n), us(n);
  std::vector<TVector3> k1x, k2x, k3x, k4x;
  std::vector<TVector3> k1u, k2u, k3u, k4u;

  GetState(x0, u0);

  // stage 1
  ComputeDerivative(k1x, k1u);
  for(i=0; i<n; i++){
    xs[i] = x0[i] + k1x[i]*(0.5*dt);
    us[i] = u0[i] + k1u[i]*(0.5*dt);
  }
  SetState(xs, us);
  // stage 2
  ComputeDerivative(k2x, k2u);
  for(i=0; i<n; i++){
    xs[i] = x0[i] + k2x[i]*(0.5*dt);
    us[i] = u0[i] + k2u[i]*(0.5*dt);
  }
  SetState(xs, us);
  // stage 3
  ComputeDerivative(k3x, k3u);
  for(i=0; i<n; i++){
    xs[i] = x0[i] + k3x[i]*dt;
    us[i] = u0[i] + k3u[i]*dt;
  }
  SetState(xs, us);
  // stage 4
  ComputeDerivative(k4x, k4u);
  for(i=0; i<n; i++){
    xs[i] = x0[i] + (k1x[i] + k2x[i]*2. + k3x[i]*2. + k4x[i]) * (dt/6.);
    us[i] = u0[i] + (k1u[i] + k2u[i]*2. + k3u[i]*2. + k4u[i]) * (dt/6.);
  }
  SetState(xs, us);
  FinishStep(dt);
}

void event::DeriveRK45(Double_t dt){
  // Dormand-Prince 5(4) (Hairer, Norsett, Wanner, Solving ODE I, II.5).
  // The 5th order solution is kept; the embedded 4th order one gives the error.
  // The last stage is evaluated at the new state, so it is the first stage of
  // the next substep (FSAL). Forces do not depend on time, so the nodes c_i are not needed.
  static const Double_t a[7][6] = {
    {0., 0., 0., 0., 0., 0.},
    {1./5., 0., 0., 0., 0., 0.},
    {3./40., 9./40., 0., 0., 0., 0.},
    {44./45., -56./15., 32./9., 0., 0., 0.},
    {19372./6561., -25360./2187., 64448./6561., -212./729., 0., 0.},
    {9017./3168., -355./33., 46732./5247., 49./176., -5103./18656., 0.},
    {35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.}};
  // 5th - 4th order weights
  static const Double_t e[7] = {71./57600., 0., -71./16695., 71./1920.,
                                -17253./339200., 22./525., -1./40.};
  // step controller : h_new = h * clamp(safe * err^(-1/5), facmin, facmax)
  static const Double_t safe = 0.9, facmin = 0.2, facmax = 5.;

  const size_t n = particles.size();
  const Double_t total = TMath::Abs(dt);
  const Double_t sgn   = (dt<0.) ? -1. : 1.;
  if(total==0. || n==0){
    FinishStep(dt);
    return;
  }
  // below this substep the error is accepted as it is (with a warning),
  // instead of looping for ever on a singular force
  const Double_t hmin = total*1e-12;
  const Bool_t   rel  = UseRelativistic();

  size_t i;
  Int_t  s, j;
  std::vector<TVector3> x0, u0, xs(n), us(n);
  std::vector<TVector3> kx[7], ku[7];

  GetState(x0, u0);
  ComputeDerivative(kx[0], ku[0]);

  Double_t h = (rk45H>0.) ? rk45H : total; // magnitude of the substep
  Double_t done = 0.;
  Bool_t   rejected = kFALSE;              // the present substep was rejected before

  Long_t   ntry = 0;                       // substeps tried in this call
  while(kTRUE){
    // A force singularity (e.g. an unsoftened attractive centre) shrinks the
    // substep without end. Inspectors and the step cap are only looked at
    // between calls, so the call itself has to give up.
    if(++ntry > RK45MAXSUBSTEPS){
      std::cerr<<"event::DeriveRK45 : more than "<<RK45MAXSUBSTEPS<<" substeps in one call of DeriveDT("
               <<dt<<"), substep "<<h<<". Integration given up (singular force ?)."<<std::endl;
      integFailed = kTRUE;
      SetState(x0, u0); // last accepted state
      break;
    }
    // the end of dt : take the rest exactly (also when only a sliver would be left)
    Double_t remain = total - done;
    Double_t htry = h;
    Bool_t   last = kFALSE;
    if(h*1.01 >= remain){
      htry = remain;
      last = kTRUE;
    }
    Double_t hs = sgn*htry;

    // stages 2..7 ; after stage 7 the particles are at the 5th order solution
    for(s=1; s<7; s++){
      for(i=0; i<n; i++){
        xs[i] = x0[i];
        us[i] = u0[i];
        for(j=0; j<s; j++){
          if(a[s][j]==0.) continue;
          xs[i] += kx[j][i] * (hs*a[s][j]);
          us[i] += ku[j][i] * (hs*a[s][j]);
        }
      }
      SetState(xs, us);
      ComputeDerivative(kx[s], ku[s]);
    }

    // RMS of the error over the moving particles, scaled by the tolerance
    Double_t sum = 0.;
    Int_t    nterm = 0;
    for(i=0; i<n; i++){
      particle* p = particles[i];
      if(p->IsInvincible()) continue;
      TVector3 ex, eu;
      for(j=0; j<7; j++){
        if(e[j]==0.) continue;
        ex += kx[j][i] * (hs*e[j]);
        eu += ku[j][i] * (hs*e[j]);
      }
      // momentum / m0 : the absolute tolerance is in c in both modes
      Double_t um = rel ? 1./p->GetM0() : 1.;
      Double_t scx = rk45Abs + rk45Rel*TMath::Max(x0[i].Mag(), xs[i].Mag());
      Double_t scu = rk45Abs + rk45Rel*TMath::Max(u0[i].Mag(), us[i].Mag())*um;
      Double_t rx = ex.Mag()/scx;
      Double_t ru = eu.Mag()*um/scu;
      sum += rx*rx + ru*ru;
      nterm += 2;
    }
    Double_t err = (nterm>0) ? TMath::Sqrt(sum/nterm) : 0.;
    if(!(err<=1e300)) err = 1e300; // NaN or inf : reject

    if(err<=1. || htry<=hmin){
      if(err>1. && !rk45Warned){
        std::cerr<<"event::DeriveRK45 : substep "<<htry<<" reached the minimum and the error is "
                 <<err<<" times the tolerance. Accepted anyway (warned once per event)."<<std::endl;
        rk45Warned = kTRUE;
      }
      rk45NAccept++;
      Double_t fac = (err>0.) ? safe*TMath::Power(err, -0.2) : facmax;
      fac = TMath::Min(facmax, TMath::Max(facmin, fac));
      if(rejected) fac = TMath::Min(fac, 1.);
      Double_t hnew = htry*fac;
      // a substep cut by the end of dt says nothing against the longer one suggested before
      if(last && htry<h && fac>=1.) hnew = TMath::Max(hnew, h);
      h = hnew;
      rejected = kFALSE;

      x0 = xs;
      u0 = us;
      kx[0] = kx[6];
      ku[0] = ku[6];
      if(last) break;
      done += htry;
    }else{
      rk45NReject++;
      h = htry*TMath::Max(facmin, safe*TMath::Power(err, -0.2));
      rejected = kTRUE;
      SetState(x0, u0);
    }
  }
  rk45H = h;
  FinishStep(dt);
}

void event::DeriveTaylor(Double_t dt){
  // Int_t pi, pj, fi;
  vector<particle*>::iterator p1_;
  vector<particle*>::iterator p2_;
  vector<force*>::iterator    ff_;
  vector<particle*>::iterator pp_;

  particle *p1, *p2;
  force *ff;

  TVector3 tempforce = TVector3();
  for(p1_=particles.begin(); p1_!=particles.end(); p1_++){
    p1 = *p1_;
    if(p1->IsInvincible()) continue;
    for(p2_=particles.begin(); p2_!=particles.end(); p2_++){
      p2 = *p2_;
      if(p1==p2){
        continue;
      }
      for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
        ff = *ff_;
        tempforce = ff->Evaluate(p1, p2);
        p1->holdForce(tempforce);
      }
    }
  }

  for(pp_=particles.begin(); pp_!=particles.end(); pp_++){
    (*pp_)->releaseForce(dt);
  }
  return;

  // Int_t i;
  // for(i=0; i<recorders.size(); i++){
  //   if(recorders[i]->IsConditioned()){
  //     recorders[i]->Save();
  //   }else{
  //     recorders[i]->JudgeSave();
  //   }
  // }
  // return;
}

void event::DeriveDTN(Int_t n, Double_t dt){
  Int_t i;
  for(i=0; i<n; i++){
    DeriveDT(dt);
    // std::cout<<"DERIVE: "<<i<<endl;
  }
}

void event::DeriveMAX(Double_t dt){
  Int_t i;
  for(i=0; i<MAXSTEPS; i++){
    DeriveDT(dt);
    // std::cout<<"DERIVE: "<<i<<endl;
  }
}

particle * event::getParticle(Int_t index){
  vector<particle*>::iterator p1_;
  particle* p1;

  Int_t i=0;
  if(index >= getNParticle()){
    std::cerr<<"CALL PARTICLE ERROR:: Index over NParticle"<<endl;
    return 0;
  }
  for(p1_=particles.begin(); p1_!=particles.end(); p1_++){
    // p1 = &p1_;
    p1 = *p1_;
    if(i==index){
      return p1;
    }
    i++;
  }
  return p1;
}

force * event::getForce(Int_t index){
  vector<force*>::iterator f1_;
  force* f1;

  Int_t i=0;
  if(index >= getNForce()){
    std::cerr<<"CALL Force ERROR:: Index over NForce"<<endl;
    return 0;
  }
  for(f1_=forces.begin(); f1_!=forces.end(); f1_++){
    f1 = *f1_;
    if(i==index){
      return f1;
    }
    i++;
  }
  return f1;
}

inspector * event::getInspector(Int_t index){
  vector<inspector*>::iterator p1_;
  inspector* p1;

  Int_t i=0;
  if(index >= getNInspector()){
    std::cerr<<"CALL INSPECTOR ERROR:: Index over NInspector"<<endl;
    return 0;
  }
  for(p1_=inspectors.begin(); p1_!=inspectors.end(); p1_++){
    // p1 = &p1_;
    p1 = *p1_;
    if(i==index){
      return p1;
    }
    i++;
  }
  return p1;
}

Long_t event::GetDeriveN(const Bool_t minima){
  
  vector<particle*>::iterator p_;
  // particle *p;

  Long_t result = minima ? std::numeric_limits<Long_t>::max() : std::numeric_limits<Long_t>::min();
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    result = minima ? TMath::Min(result, (*p_)->GetECNT()) : TMath::Max(result, (*p_)->GetECNT());
  }
  return result;

}

Double_t event::GetNetKE(const Bool_t SR){

  vector<particle*>::iterator p_;
  particle *p;

  Double_t tempKE = 0.;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    p = *p_;
    if(p->IsInvincible()) continue;
    tempKE += (0.5)*(p->GetM(!SR))*(p->GetV().Mag2());
  }
  return tempKE;
}

Double_t event::GetNetKERel(){
  vector<particle*>::iterator p_;
  Double_t tempKE = 0.;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    particle* p = *p_;
    if(p->IsInvincible()) continue;
    // gamma - 1 = v^2 / (s (1 + s)), s = sqrt(1 - v^2) : no cancellation at low speed
    Double_t v2 = p->GetV().Mag2();
    Double_t s  = TMath::Sqrt(1. - v2);
    tempKE += p->GetM0() * v2 / (s*(1. + s));
  }
  return tempKE;
}

TVector3 event::GetNetMomentum(){
  const Bool_t rel = UseRelativistic();
  vector<particle*>::iterator p_;
  TVector3 result = TVector3();
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    particle* p = *p_;
    if(p->IsInvincible()) continue;
    TVector3 v = p->GetV();
    result += rel ? v * (p->GetM0()/TMath::Sqrt(1. - v.Mag2())) : v * p->GetM0();
  }
  return result;
}

TVector3 event::GetNetAngularMomentum(TVector3 origin){
  const Bool_t rel = UseRelativistic();
  vector<particle*>::iterator p_;
  TVector3 result = TVector3();
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    particle* p = *p_;
    if(p->IsInvincible()) continue;
    TVector3 v = p->GetV();
    TVector3 mom = rel ? v * (p->GetM0()/TMath::Sqrt(1. - v.Mag2())) : v * p->GetM0();
    result += (p->GetX() - origin).Cross(mom);
  }
  return result;
}

Double_t event::GetPE(Int_t index, Bool_t OnlyInvincibles){
  return GetPE(getParticle(index), OnlyInvincibles);
}

Double_t event::GetPE(particle* p, Bool_t OnlyInvincibles){
  vector<particle*>::iterator pt_;
  vector<force*>::iterator    ff_;

  particle *pt;
  force *ff;

  Double_t tempPE = 0.;

  for(pt_=particles.begin(); pt_!=particles.end(); pt_++){
    pt = *pt_;
    if(pt==p){
      continue;
    }
    for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
      ff = *ff_;
      if((!(pt->IsInvincible()))&&OnlyInvincibles)continue;
      tempPE += ff->Potential(p, pt);
    }
  }
  return tempPE;
}

Double_t event::GetNetPE(){

  // vector<particle*>::iterator p1_;
  // vector<particle*>::iterator p2_;
  // vector<force*>::iterator    ff_;

  // particle *p1, *p2;
  // force *ff;

  // Double_t tempPE = 0.;
  // for(p1_=particles.begin(); p1_!=particles.end(); p1_++){
  //   p1 = *p1_;
  //   for(p2_=particles.begin(); p2_!=particles.end(); p2_++){
  //     p2 = *p2_;
  //     if(p1==p2){
  //       continue;
  //     }
  //     for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
  //       ff = *ff_;
        
  //       tempPE += ff->Potential(p1, p2);
  //     }
  //   }
  // }
  // return tempPE;

  Int_t i;
  Double_t tempPE = 0.;
  for(i=0; i<getNParticle(); i++){
    // std::cout<<"GETPARTICLEPE "<<i<<std::endl;
    tempPE += 0.5*(GetPE(i));
  }
  return tempPE;
}


void event::preDerive(){
  if(preDerived){
    delete pEvent;
    pEvent = nullptr;
    // cout<<"Clean previous pre-derive"<<endl;
  }
  pEvent = new event(*this);
  preDerived = kTRUE;
  pEvent->DeriveInspect(1);
  // cout<<"NEW EVENT COMP"<<endl;
  return;
}

Bool_t event::resetPreDerive(){
  if(preDerived){
    delete pEvent;
    pEvent = nullptr;
    preDerived=kFALSE;
    return kTRUE;
  }else{
    return kFALSE;
  }
}

Bool_t event::DeriveInspect(Int_t iperiod, Double_t dt){
  lastStop = kStopNone;
  if(iperiod<1) iperiod = 1;

  // An inspector with an undefined method never gets satisfied : refuse to start.
  vector<inspector*>::iterator i_;
  for(i_=inspectors.begin(); i_!=inspectors.end(); i_++){
    if(!(*i_)->IsValid()){
      std::cerr<<"event::DeriveInspect : inspector with undefined method / "
               <<(*i_)->GetMethod()<<" . Nothing derived."<<std::endl;
      lastStop = kStopInvalid;
      return kFALSE;
    }
  }

  Long_t i;
  for(i=0; ; i++){
    if(i%iperiod==0 && Inspect()){
      if(lastStop!=kStopAbort) lastStop = kStopInspector;
      break;
    }
    if(maxInspectSteps>0 && i>=maxInspectSteps){
      std::cerr<<"event::DeriveInspect : stopped by the maximum number of steps ("
               <<maxInspectSteps<<"), inspectors are not satisfied."<<std::endl;
      lastStop = kStopMaxSteps;
      break;
    }
    DeriveDT(dt);
    if(integFailed){
      lastStop = kStopFailed;
      break;
    }
  }
  return lastStop==kStopInspector;
}

Bool_t event::Inspect(){
  // 1. Abort inspectors : any of them stops the event, whatever the others say.
  // 2. The other inspectors are combined by the inspect mode (AND by default).
  vector<inspector*>::iterator i_;
  inspector* i;

  Int_t  nnormal = 0, nabort = 0;
  Bool_t resultAND = kTRUE;
  Bool_t resultOR  = kFALSE;

  for(i_=inspectors.begin(); i_!=inspectors.end(); i_++){
    i = *i_;
    if(i->IsAbort()){
      nabort++;
      if(i->Inspect()){
        lastStop = kStopAbort;
        return kTRUE;
      }
      continue;
    }
    nnormal++;
    Bool_t r = i->Inspect();
    resultAND &= r;
    resultOR  |= r;
  }

  if(nnormal==0){
    // only abort inspectors : keep going until one of them fires.
    // no inspector at all : nothing to wait for (as before).
    return nabort==0;
  }
  return imode==kInspectOR ? resultOR : resultAND;
}

Double_t event::GetTime(){
  vector<particle*>::iterator p_;
  Double_t result = 0.;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    result = TMath::Max(result, (*p_)->GetTime());
  }
  return result;
}

Double_t event::GetMinDistance(){
  Double_t result = std::numeric_limits<Double_t>::max();
  size_t i, j;
  for(i=0; i<particles.size(); i++){
    for(j=i+1; j<particles.size(); j++){
      result = TMath::Min(result, (particles[i]->GetX()-particles[j]->GetX()).Mag());
    }
  }
  return result;
}

Double_t event::GetMaxAcceleration(){
  std::vector<TVector3> acc;
  ComputeAcceleration(acc);
  Double_t result = 0.;
  size_t i;
  for(i=0; i<acc.size(); i++){
    result = TMath::Max(result, acc[i].Mag());
  }
  return result;
}

Bool_t event::IsInsideVolume(Bool_t all){
  // Without a volume, everything is regarded as inside.
  if(!volume) return kTRUE;
  vector<particle*>::iterator p_;
  Bool_t resultAll = kTRUE;
  Bool_t resultAny = kFALSE;
  Int_t  n = 0;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    if((*p_)->IsInvincible()) continue;
    n++;
    Bool_t in = volume->Inside((*p_)->GetX());
    resultAll &= in;
    resultAny |= in;
  }
  if(n==0) return kTRUE;
  return all ? resultAll : resultAny;
}

void event::preDeriveSetup(){
  // this->getInspector(0)->GetParticle1()->



  // return;
}

TString event::Print(Bool_t onlymechanic, Bool_t mute, Bool_t pprint){
  TString result = TString("[ Event ");
  
  if(pprint){
    result = result + "\n";
  }

  if(particles.size()>0){
    result += TString("| Particles ");
    result += TString("{ ");
    std::vector<particle*>::iterator i;
    for(i=particles.begin(); i<particles.end(); i++){
      result += (*i)->Print(onlymechanic, kTRUE);
    }
    result += TString(" }");
    if(pprint){
      result = result + "\n";
    }
  }

  if(inspectors.size()>0){
    result += TString("| Inspectors ");
    result += TString("{ ");
    if(pprint){
      result += "\n";
    }
    std::vector<inspector*>::iterator i;
    for(i=inspectors.begin(); i<inspectors.end(); i++){
      result += (*i)->Print(kTRUE, kTRUE, kTRUE, pprint);
    }
    result += TString(" }");
  }

  result += TString(" ]");

  if(!mute){
    cout<<result<<endl; 
  }

  return result;
}