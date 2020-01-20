#ifndef H_EXPERIMENT
#define H_EXPERIMENT 1;

#include "TROOT.h"
#include "particle.hpp"
#include "event.hpp"

#include "iostream"
#include "cstdarg"

class experiment{
public:
  experiment(){}
  experiment(event* eventT_);
  void SetTemplate(event* eventT_){eventT=eventT_;}
  void MakeEvent(Int_t ninparticle_, ...);
  void OffEvent(){eventOn=false;}
  void DelEvent(){eventOn=false; delete currentEvent;}
  event* GetEvent(){return currentEvent;}

  friend class experimentP;
protected:
private:
  event * eventT;
  Bool_t eventOn = false;
  event * currentEvent;
};

class experimentP: public experiment{
public:
  experimentP(Int_t ninparticle_){ninparticle=ninparticle_;}
  void SetNInParticle(Int_t ninparticle_){ninparticle = ninparticle_;}
  Int_t GetNInParticle(){return ninparticle;}
  void MakeEvent(particle* p1, ...);
private:
  Int_t ninparticle=1;
};


#endif