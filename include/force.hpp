#ifndef H_FORCE
#define H_FORCE 1;

#include "TROOT.h"
#include "TVector3.h"
#include "particle.hpp"

class force{
public:
  virtual TVector3 Evaluate(particle* p1_, particle* p2_){
    TVector3 a_;
    std::cout<<"NULL Force"<<endl; return a_;
  }
  virtual Double_t Potential(particle* p1_, particle* p2_){
    Double_t a_ = 0.;
    std::cout<<"NULL Force, Potential"<<endl; return a_;
  };
protected:
  force(){}
private:
};

class coulombForce: public force{
public:
  coulombForce():force(){}
  TVector3 Evaluate(particle* p1_, particle* p2_);
  Double_t Potential(particle* p1_, particle* p2_);
protected:
Double_t coulombK = 0.00154629;
private:
};

//G = Quantity[1.2331*10^-17, ("Femtometers" ("SpeedOfLight")^2)/("AtomicMassUnit")]
class universalForce: public force{ //Gravity
public:
  universalForce():force(){}
  TVector3 Evaluate(particle* p1_, particle* p2_);
  Double_t Potential(particle* p1_, particle* p2_);
protected:
private:
};

#endif
