#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TVector3.h"

#include "particle.hpp"
#include "force.hpp"

using namespace std;

TVector3 coulombForce::Evaluate(particle* p1, particle* p2){
  // Double_t coulombK = 0.00154629;
  Double_t len, q1, q2;
  q1 = p1->GetQ();
  q2 = p2->GetQ();
  if(softening!=0.){
    // F = k q1 q2 d / (r^2 + eps^2)^(3/2) = -grad of the softened potential
    TVector3 d = p1->GetX() - p2->GetX();
    Double_t s2 = d.Mag2() + softening*softening;
    return TVector3(d * (coulombK * q1 * q2 / (s2*sqrt(s2))));
  }
  // Unsoftened : kept exactly as before so that default results do not change.
  len = (p1->GetX() - p2->GetX()).Mag();
  TVector3 unit = (p1->GetX() - p2->GetX())*(1/len);
  TVector3 tempforce = TVector3(unit * coulombK * q1 * q2 * (1/pow(len,2)));
  return tempforce;
}

Double_t coulombForce::Potential(particle* p1, particle* p2){
  Double_t len, q1, q2;
  q1 = p1->GetQ();
  q2 = p2->GetQ();
  if(softening!=0.){
    Double_t s2 = ((p1->GetX())-(p2->GetX())).Mag2() + softening*softening;
    return coulombK * q1 * q2 / sqrt(s2);
  }
  len = ((p1->GetX())-(p2->GetX())).Mag();
  Double_t temppotential = coulombK * q1 * q2 * (1/len);
  return temppotential;
}

TVector3 universalForce::Evaluate(particle* p1, particle* p2){
  // Rest masses : the integrators are non-relativistic (event uses GetM(kTRUE)).
  Double_t m1 = p1->GetM0();
  Double_t m2 = p2->GetM0();
  TVector3 d = p1->GetX() - p2->GetX();
  Double_t s2 = d.Mag2() + softening*softening;
  // F = -G m1 m2 d / s^3 : towards p2
  return TVector3(d * (-gravityG * m1 * m2 / (s2*sqrt(s2))));
}

Double_t universalForce::Potential(particle* p1, particle* p2){
  Double_t m1 = p1->GetM0();
  Double_t m2 = p2->GetM0();
  Double_t s2 = ((p1->GetX())-(p2->GetX())).Mag2() + softening*softening;
  return -gravityG * m1 * m2 / sqrt(s2);
}

TVector3 yukawaForce::Evaluate(particle* p1, particle* p2){
  // U = -g e^{-r/a} / r
  // dU/dr = g e^{-r/a} (1/(a r) + 1/r^2), F = -dU/dr * d/r
  TVector3 d = p1->GetX() - p2->GetX();
  Double_t r = d.Mag();
  Double_t dUdr = g * exp(-r/a) * (1./(a*r) + 1./(r*r));
  return TVector3(d * (-dUdr / r));
}

Double_t yukawaForce::Potential(particle* p1, particle* p2){
  Double_t r = ((p1->GetX())-(p2->GetX())).Mag();
  return -g * exp(-r/a) / r;
}

// f(1-f) for the Fermi function f = 1/(1+e^x), without overflow of e^x
// on either side (f(1-f) is even in x).
static Double_t fermiDerivative(Double_t x){
  Double_t e = exp(-fabs(x));
  return e / ((1.+e)*(1.+e));
}

TVector3 woodsSaxonForce::Evaluate(particle* p1, particle* p2){
  // U = -V0 f, f = 1/(1+e^x), x = (r-R)/a
  // dU/dr = V0 f(1-f) / a, F = -dU/dr * d/r
  TVector3 d = p1->GetX() - p2->GetX();
  Double_t r = d.Mag();
  if(r==0.) return TVector3(); // direction undefined, magnitude finite
  Double_t dUdr = V0 * fermiDerivative((r-R)/a) / a;
  return TVector3(d * (-dUdr / r));
}

Double_t woodsSaxonForce::Potential(particle* p1, particle* p2){
  Double_t r = ((p1->GetX())-(p2->GetX())).Mag();
  Double_t x = (r-R)/a;
  // 1/(1+e^x) written with e^{-x} for large x, so that e^x does not overflow
  if(x>0.){
    Double_t e = exp(-x);
    return -V0 * e / (1.+e);
  }
  return -V0 / (1.+exp(x));
}
