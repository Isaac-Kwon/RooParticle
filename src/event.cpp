#include "iostream"
#include "list"
#include "stdio.h"

#include "TROOT.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "TVectorT.h"

#include "force.hpp"

#include "particle.hpp"
#include "event.hpp"

event::event(){
}

void event::DeriveDT(Double_t dt){
  // Int_t pi, pj, fi;
  list<particle*>::iterator p1_;
  list<particle*>::iterator p2_;
  list<force*>::iterator    ff_;

  list<particle*>::iterator pp_;
  particle* pp;

  particle* p1;
  particle* p2;
  force* ff;

  TVectorD tempforce = TVectorD(3);

  for(p1_=particles.begin(); p1_!=particles.end(); p1_++){
    // p1 = &p1_;
    p1 = *p1_;
    for(p2_=particles.begin(); p2_!=particles.end(); p2_++){
      // p2 = &p2_;
      p2 = *p2_;
      if(p1==p2){
        continue;
      }
      for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
        // ff = &ff_;
        // std::cout<<"FORCEPSTART"<<endl;
        ff = *ff_;
        // std::cout<<"FORCESTART"<<endl;
        // ff->Evaluate(p1, p2).Print();
        tempforce = ff->Evaluate(p1, p2);
        // std::cout<<"HOLD FORCE"<<endl;
        // tempforce.Print();
        p1->holdForce(tempforce);
      }
    }
  }

  for(pp_=particles.begin(); pp_!=particles.end(); pp_++){
    pp = *pp_;
    pp->releaseForce(dt);
  }
}
