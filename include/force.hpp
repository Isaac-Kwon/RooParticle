#ifndef H_FORCE
#define H_FORCE 1;

#include "TVectorT.h"
#include "TVectorD.h"
#include "particle.hpp"

using namespace std;

class force{
public:
  virtual TVectorD Evaluate(particle* p1_, particle* p2_){
    TVectorD a_ = TVectorD(3); cout<<"NULL FORCE"<<endl; return a_;
  }
protected:
  force(){}
private:
};

class coulombForce: public force{
public:
  coulombForce():force(){}
  TVectorD Evaluate(particle* p1_, particle* p2_);
protected:
private:
};

#endif
