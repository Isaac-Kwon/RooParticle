#ifndef H_EVENT
#define H_EVENT 1;

#include "iostream"
#include "stdio.h"
#include "list"

#include "TROOT.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "TVectorT.h"

#include "particle.hpp"
#include "force.hpp"
#include "eventVolume.hpp"

using namespace std;

//class for one event.
class event{
public:
  event();

  void DeriveDT(Double_t dt=1);
  void AddParticle(particle* particle_){
    particles.push_back(particle_); nparticle++;
  }
  void AddForce(force * force_){forces.push_back(force_); nforce++;}
  void SetVolume(eventVolume * volume_){volume = volume_;}

  Int_t getNParticle();
  Int_t getNForce();

protected:
private:
  list<particle*> particles;
  list<force*> forces;
  eventVolume * volume;

  Int_t nparticle = 0;
  Int_t nforce = 0;
};

#endif
