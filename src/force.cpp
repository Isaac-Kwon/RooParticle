#include "iostream"
#include <math.h>

#include "TROOT.h"
#include "TVectorT.h"

#include "particle.hpp"
#include "force.hpp"

using namespace std;

TVectorD coulombForce::Evaluate(particle* p1_, particle* p2_){
  EMparticle* p1 = (EMparticle *) p1_;
  EMparticle* p2 = (EMparticle *) p2_;
  Double_t coulombK = 0.00154629;
  Double_t len, q1, q2;
  len = sqrt((p1->GetX() - p2->GetX()).Norm2Sqr());
  TVectorD unit = (p1->GetX() - p2->GetX())*(1/len);
  q1 = p1->GetQ();
  q2 = p2->GetQ();
  //  = TVectorD(3);
  TVectorD tempforce = TVectorD(unit * coulombK * q1 * q2 * (1/pow(len,2)));
  // tempforce.Print();
  // std::cout<<"QQ"<<endl;
  return tempforce;
}
