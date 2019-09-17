#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TVector3.h"

#include "particle.hpp"
#include "force.hpp"

using namespace std;

TVector3 coulombForce::Evaluate(particle* p1_, particle* p2_){
  EMparticle* p1 = (EMparticle *) p1_;
  EMparticle* p2 = (EMparticle *) p2_;
  Double_t coulombK = 0.00154629;
  Double_t len, q1, q2;
  len = (p1->GetX() - p2->GetX()).Mag();
  TVector3 unit = (p1->GetX() - p2->GetX())*(1/len);
  q1 = p1->GetQ();
  q2 = p2->GetQ();
  TVector3 tempforce = TVector3(unit * coulombK * q1 * q2 * (1/pow(len,2)));
  return tempforce;
}
