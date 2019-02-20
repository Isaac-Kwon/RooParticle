#include "iostream"
#include "list"
#include "stdio.h"

#include "RPConfig.hpp"

#include "TROOT.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "TVectorT.h"

#include "force.hpp"

#include "particle.hpp"
#include "event.hpp"

//General Constructor (Not Include any procedure)
event::event(){
}

//Copy Constructor
event::event(event &e){
  list<particle*>::iterator p_;
  particle* p;
  particle* p2;
  for(p_=e.particles.begin(); p_!=e.particles.end(); p_++){
    p = *p_;
    p2 = new particle(*p);
    AddParticle(p2);
  }

  //copy particles
  list<force*>::iterator ff_;
  force* ff;
  for(ff_=e.forces.begin(); ff_!=e.forces.end(); ff_++){
    ff = *ff_;
    AddForce(ff);
  }

  //copy forces
  volume = e.volume;

  nparticle = e.nparticle;
  nforce = e.nforce;
}

event::~event(){

  list<particle*>::iterator p_;
  particle* p;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    p = *p_;
    delete p;
  }

  //copy particles
  // list<force*>::iterator ff_;
  // force* ff;
  // for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
  //   ff = *ff_;
  //   delete ff;
  // }

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

void event::DeriveMAX(Double_t dt){
  Int_t i;
  for(i=0; i<MAXSTEPS; i++){
    DeriveDT(dt);
    // std::cout<<"DERIVE: "<<i<<endl;
  }
}
