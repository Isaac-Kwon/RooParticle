
#ifndef H_FORCE
#define H_FORCE 1;

#include "TVectorT.h"
#include "particle.hpp"

class force{
public:
  virtual TVectorD Evaluate(particle* p1, particle* p2){TVectorD a_ = TVectorD(3); return a_;}
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
