#include "iostream"

#include "TROOT.h"
#include "TVectorT.h"
#include "TVector3.h"

#include "particle.hpp"

#include <math.h>

using namespace std;

// particle::particle(){
// }

particle::particle(Double_t m0_,
                   // Double_t q_  = 1,
                   TVector3 x0_,
                   TVector3 v0_,
                   Bool_t   Invincible_,
                   Bool_t   Record_){
 m0 = m0_;
 // q  = q_;
 x  = x0_;
 v  = v0_;
 hF = TVector3();
 Invincible = Invincible_;
 if(Record_){
   ppath = new path();
   Record = true;
 }else{
   Record = false;
 }
}

//copy constructor
particle::particle(const particle & p){
  name       = p.name;
  x          = p.x;
  v          = p.v;
  m0         = p.m0;
  Invincible = p.Invincible;
  Record     = p.Record;
  hold       = p.hold;
  hF         = p.hF;
  etime      = p.etime;
  q          = p.q;
  if(p.Record){
   ppath = new path(*(p.ppath));
   Record = kTRUE;
  }else{
   Record = kFALSE;
 }
}

void particle::holdForce(const TVector3 f){
  if(!hold){
    hold = true;
    hF = TVector3();
    // std::cout<<"ZERO"<<endl;
    // hF.Print();
  }
  hF+=f;
}

void particle::releaseForce(const Double_t dt){
  applyForce(hF, dt);
  resetHold();
}

void particle::resetHold(){
  hold=false;
  hF = TVector3();
}

void particle::applyForce(const TVector3 f, const Double_t dt){
  if(Invincible==true){
    if(Record==true){
      ppath->AddRecord(x,v, etime);
    }
    return;
  }
  Double_t m_ = GetM(false);
  TVector3 a = f * (1/m_);
  applyDX(
    a * (0.5) * pow(dt,2) + v * dt
  );
  applyDV(
    a * dt
  );
  etime += dt;
  if(Record==true){
    // cout<<"RECORD"<<endl;
    ppath->AddRecord(x,v, etime);
    // cout<<"RECORD E"<<endl;
    }
  // cout<<"APPLYFORCE"<<endl;
  // f.Print();
  // cout<<dt<<" / "<<etime<<endl;
}

Double_t particle::GetM(Bool_t restmass){
  if(restmass){
  return m0;
  }
  Double_t v_ = v.Mag();
  return ( 1/sqrt( 1-pow(v_,2) )) * m0;
}

TString onelineTVectorS(TVector3 vec){
  TString result = TString("");
  result += TString("[");
  Int_t i;
  for(i=0; i<2; ++i){
    // cout<<i<<endl;
    result+=vec(i);
    result+=", ";
  }
  result += vec(i);
  result += TString("]");
  return result;
}

TString particle::Print(Bool_t onlymechanic, Bool_t mute){
  TString result = TString("");
  result += TString("[ Particle ");
  if(!onlymechanic){
  }
  result += TString("x: ");
  result += onelineTVectorS(GetX());
  result += TString(", ");
  result += TString("v: ");
  result += onelineTVectorS(GetV());
  result += TString(" ]");

  if(!mute){
    cout<<result<<endl;
    

  }

  return result;

  
}