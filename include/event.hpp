#ifndef H_EVENT
#define H_EVENT 1;

#include "iostream"
#include "stdio.h"
#include "list"

#include "TROOT.h"
#include "TString.h"
#include "TMath.h"


#include "particle.hpp"
#include "force.hpp"
#include "eventVolume.hpp"
#include "inspector.hpp"
#include "recorder.hpp"

class inspector;

//class for one event.
class event{
public:
  event();
  event(event &e);
  virtual ~event();

  //Setup Methods
  particle* AddParticle(particle* particle_){particles.push_back(particle_); nparticle++; return particle_;}
  void AddForce(force * force_){forces.push_back(force_); nforce++;}
  void AddInspector(inspector* inspector_){inspectors.push_back(inspector_);}
  void SetVolume(eventVolume * volume_){volume = volume_;}

  //DerivingMethod
  virtual void DeriveDT(Double_t dt = 1);
  void DeriveDTN(Int_t n = 1, Double_t dt = 1);
  void DeriveMAX(Double_t dt = 1);
  void DeriveInspect(Int_t iperiod=1);

  //GetMethod
  particle * getParticle(Int_t index);
  force * getForce(Int_t index);
  inspector * getInspector(Int_t index);
  Int_t getNParticle(){return nparticle;}
  Int_t getNForce(){return nforce;}
  Int_t getNInspector(){return inspectors.size();}

  //Calculation
  Double_t GetNetKE();
  Double_t GetNetPE();
  Double_t GetNetE(){return GetNetKE()+GetNetPE();}
  

  //Inspection Deriving and Pre-Deriving.
  void preDerive();
  void preDeriveSetup();
  Bool_t resetPreDerive();
  Bool_t Inspect();
  TString Print(Bool_t onlymechanic=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);
  Bool_t SetInitial();

protected:
  std::vector<particle*> particles;
  std::vector<force*> forces;

  Int_t nparticle = 0;
  Int_t nforce = 0;

private:
  eventVolume * volume;
  event * pEvent;
  Bool_t preDerived = kFALSE;
  
  std::vector<inspector*> inspectors;
};

#endif
