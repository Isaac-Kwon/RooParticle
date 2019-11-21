#include "iostream"
#include "list"
#include "vector"
#include "stdio.h"

#include "RPConfig.hpp"

#include "TROOT.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "TVector3.h"

#include "force.hpp"

#include "particle.hpp"
#include "event.hpp"
#include "inspector.hpp"

event::event(){
}

//Copy Constructor
event::event(event &e){
  //Copy all particle into new event
  vector<particle*>::iterator p_;
  for(p_=e.particles.begin(); p_!=e.particles.end(); p_++){
    AddParticle(new particle(**p_));
  }

  //Add all force into new event (not duplicate)
  vector<force*>::iterator ff_;
  for(ff_=e.forces.begin(); ff_!=e.forces.end(); ff_++){
    AddForce(*ff_);
  }

  //copy inspectors : link all particles in each inspectors
  vector<inspector*>::iterator ii_;

  Int_t ip1=0, ip2=0;
  for(ii_=e.inspectors.begin(); ii_!=e.inspectors.end(); ii_++){
    if((*ii_)->IsInspectorP()){
      inspectorP* iis; //source inspector
      iis = (inspectorP*) *ii_;
      Int_t i;
      for(i=0; i<e.particles.size();i++){
        if((iis->GetParticle1())==e.particles[i]){
          ip1 = i;
        }
        if((iis->GetParticle2())==e.particles[i]){
          ip2 = i;
        }
      }
      AddInspector(new inspectorP( // Make new inspector (linking particles)
        getParticle(ip1), getParticle(ip2),
        iis->GetMethod(), iis->GetCriterion())
        );
    }else if((*ii_)->IsInspectorE()){
      inspectorE* iie = (inspectorE*) *ii_;
      AddInspector(new inspectorE(iie->GetEvent(), iie->GetMethod(), iie->GetCriterion()));
    }
  }

  //copy other properties
  volume = e.volume;
  nparticle = e.nparticle;
  nforce = e.nforce;
}

event::~event(){
  // Delete all particle
  // Inspectors are automatically deleted.
  vector<particle*>::iterator p_;
  particle* p;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    p = *p_;
    delete p;
  }
}

void event::DeriveDT(Double_t dt){
  // Int_t pi, pj, fi;
  vector<particle*>::iterator p1_;
  vector<particle*>::iterator p2_;
  vector<force*>::iterator    ff_;
  vector<particle*>::iterator pp_;

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
  }

  for(pp_=particles.begin(); pp_!=particles.end(); pp_++){
    (*pp_)->releaseForce(dt);
  }
  return;

  // Int_t i;
  // for(i=0; i<recorders.size(); i++){
  //   if(recorders[i]->IsConditioned()){
  //     recorders[i]->Save();
  //   }else{
  //     recorders[i]->JudgeSave();
  //   }
  // }
  // return;
}

void event::DeriveDTN(Int_t n, Double_t dt){
  Int_t i;
  for(i=0; i<n; i++){
    DeriveDT(dt);
    // std::cout<<"DERIVE: "<<i<<endl;
  }
}

void event::DeriveMAX(Double_t dt){
  Int_t i;
  for(i=0; i<MAXSTEPS; i++){
    DeriveDT(dt);
    // std::cout<<"DERIVE: "<<i<<endl;
  }
}

particle * event::getParticle(Int_t index){
  vector<particle*>::iterator p1_;
  particle* p1;

  Int_t i=0;
  if(index >= getNParticle()){
    std::cerr<<"CALL PARTICLE ERROR:: Index over NParticle"<<endl;
    return 0;
  }
  for(p1_=particles.begin(); p1_!=particles.end(); p1_++){
    // p1 = &p1_;
    p1 = *p1_;
    if(i==index){
      return p1;
    }
    i++;
  }
  return p1;
}

force * event::getForce(Int_t index){
  vector<force*>::iterator f1_;
  force* f1;

  Int_t i=0;
  if(index >= getNForce()){
    std::cerr<<"CALL Force ERROR:: Index over NForce"<<endl;
    return 0;
  }
  for(f1_=forces.begin(); f1_!=forces.end(); f1_++){
    f1 = *f1_;
    if(i==index){
      return f1;
    }
    i++;
  }
  return f1;
}

inspector * event::getInspector(Int_t index){
  vector<inspector*>::iterator p1_;
  inspector* p1;

  Int_t i=0;
  if(index >= getNInspector()){
    std::cerr<<"CALL INSPECTOR ERROR:: Index over NInspector"<<endl;
    return 0;
  }
  for(p1_=inspectors.begin(); p1_!=inspectors.end(); p1_++){
    // p1 = &p1_;
    p1 = *p1_;
    if(i==index){
      return p1;
    }
    i++;
  }
  return p1;
}

Double_t event::GetNetKE(const Bool_t SR){

  vector<particle*>::iterator p_;
  particle *p;

  Double_t tempKE = 0.;
  for(p_=particles.begin(); p_!=particles.end(); p_++){
    p = *p_;
    if(p->IsInvincible()) continue;
    tempKE += (0.5)*(p->GetM(!SR))*(p->GetV().Mag2());
  }
  return tempKE;
}

Double_t event::GetPE(Int_t index, Bool_t OnlyInvincibles){
  return GetPE(getParticle(index), OnlyInvincibles);
}

Double_t event::GetPE(particle* p, Bool_t OnlyInvincibles){
  vector<particle*>::iterator pt_;
  vector<force*>::iterator    ff_;

  particle *pt;
  force *ff;

  Double_t tempPE = 0.;

  for(pt_=particles.begin(); pt_!=particles.end(); pt_++){
    pt = *pt_;
    if(pt==p){
      continue;
    }
    for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
      ff = *ff_;
      if((!(pt->IsInvincible()))&&OnlyInvincibles)continue;
      tempPE += ff->Potential(p, pt);
    }
  }
  return tempPE;
}

Double_t event::GetNetPE(){

  // vector<particle*>::iterator p1_;
  // vector<particle*>::iterator p2_;
  // vector<force*>::iterator    ff_;

  // particle *p1, *p2;
  // force *ff;

  // Double_t tempPE = 0.;
  // for(p1_=particles.begin(); p1_!=particles.end(); p1_++){
  //   p1 = *p1_;
  //   for(p2_=particles.begin(); p2_!=particles.end(); p2_++){
  //     p2 = *p2_;
  //     if(p1==p2){
  //       continue;
  //     }
  //     for(ff_=forces.begin(); ff_!=forces.end(); ff_++){
  //       ff = *ff_;
        
  //       tempPE += ff->Potential(p1, p2);
  //     }
  //   }
  // }
  // return tempPE;

  Int_t i;
  Double_t tempPE = 0.;
  for(i=0; i<getNParticle(); i++){
    tempPE += GetPE(i);
  }
  return tempPE;
}


void event::preDerive(){
  if(preDerived){
    delete pEvent;
    // cout<<"Clean previous pre-derive"<<endl;
  }
  pEvent = new event(*this);
  preDerived = kTRUE;
  pEvent->DeriveInspect(1);
  // cout<<"NEW EVENT COMP"<<endl;
  return;
}

Bool_t event::resetPreDerive(){
  if(preDerived){
    delete pEvent;
    preDerived=kFALSE;
    return kTRUE;
  }else{
    return kFALSE;
  }
}

void event::DeriveInspect(Int_t iperiod){
  Int_t i;
  for(i=0; i%iperiod==0 ? !Inspect() : kTRUE ;i++){
    DeriveDT();
  }
  return;
}

Bool_t event::Inspect(){
  vector<inspector*>::iterator i_;
  inspector* i;
  Bool_t result=kTRUE;

  for(i_=inspectors.begin(); i_!=inspectors.end(); i_++){
    i = *i_;
    result &= i->Inspect();
  }
  return result; 
}

void event::preDeriveSetup(){
  // this->getInspector(0)->GetParticle1()->



  // return;
}

TString event::Print(Bool_t onlymechanic, Bool_t mute, Bool_t pprint){
  TString result = TString("[ Event ");
  
  if(pprint){
    result = result + "\n";
  }

  if(particles.size()>0){
    result += TString("| Particles ");
    result += TString("{ ");
    std::vector<particle*>::iterator i;
    for(i=particles.begin(); i<particles.end(); i++){
      result += (*i)->Print(onlymechanic, kTRUE);
    }
    result += TString(" }");
    if(pprint){
      result = result + "\n";
    }
  }

  if(inspectors.size()>0){
    result += TString("| Inspectors ");
    result += TString("{ ");
    if(pprint){
      result += "\n";
    }
    std::vector<inspector*>::iterator i;
    for(i=inspectors.begin(); i<inspectors.end(); i++){
      result += (*i)->Print(kTRUE, kTRUE, kTRUE, pprint);
    }
    result += TString(" }");
  }

  result += TString(" ]");

  if(!mute){
    cout<<result<<endl; 
  }

  return result;
}