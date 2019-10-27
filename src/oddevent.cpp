#include "TROOT.h"
#include "oddevent.hpp"
#include "event.hpp"

void oddevent::DeriveDT(Double_t dt){
  // Int_t pi, pj, fi;
  vector<particle*>::iterator p1_;
  vector<particle*>::iterator p2_;
  vector<force*>::iterator    ff_;

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
    p1->releaseForce(dt);
  }
  return;
}