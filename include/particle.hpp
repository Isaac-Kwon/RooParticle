#ifndef H_PARTICLE
#define H_PARTICLE 1;

#include "TROOT.h"
#include "TString.h"
#include "TVector3.h"

#include "path.hpp"

class inspector;

// using namespace std;

class particle{
public:
  // particle();
  particle(Double_t m0_ = 1,
           //Double_t q_  = 1,
           TVector3 x0_ = TVector3(),
           TVector3 v0_ = TVector3(),
           Bool_t   Invincible_ = false,
           Bool_t   Record_     = false);
  //Force Methods (Holding and Release)
  particle(const particle & p); //Copy Constructor
  ~particle(){if(Record)delete ppath;} //Destructor (NO Pointer Variable in class)

  void holdForce(const TVector3 f);
  void releaseForce(const Double_t dt);
  void resetHold();

  //Applying Methods
  void applyMechanic(const TVector3 x_, const TVector3 v_){x = x_; v = v_;}
  void applyX(const TVector3 x_){x=x_;}
  void applyV(const TVector3 v_){v=v_;}
  void applyForce(const TVector3 f, const Double_t dt=1.);
  void applyDX(const TVector3 dx){x += dx;}
  void applyDV(const TVector3 dv){v += dv;}
  void applyDelta(const TVector3 dx, const TVector3 dv){x += dx; v += dv;}

  void SetName(const TString name_){name=name_;}
  void AddRecord(){if(Record)ppath->AddRecord(x,v,etime);}

  TVector3 GetX(){return TVector3(x);}
  TVector3 GetV(){return TVector3(v);}
  TVector3 GetP(){return TVector3(v * GetM());}
  Long_t   GetECNT(){return ecount;} //Eveoled number of point
  Double_t GetM0(){return m0;}
  Double_t GetM(Bool_t restmass=false);
  Double_t GetTime(){return etime;}
  path*    GetPath(){return ppath;}

  Bool_t   IsInvincible(){return Invincible;}
  Bool_t   IsRecording(){return Record;}

  TString  Print(Bool_t onlymechanic=kFALSE, Bool_t mute=kFALSE);

  friend class inspector;

protected:
  Double_t q          = 0;
  
private:
  TString  name;
  TVector3 x          = TVector3();
  TVector3 v          = TVector3();
  Double_t m0         = 1;

  Bool_t   Invincible = false;
  Bool_t   hold       = false;
  TVector3 hF         = TVector3();

  Double_t etime      = 0.;
  Long_t   ecount      = 0;

  Bool_t   Record     = false;
  path*    ppath;
};

class EMparticle : public particle{
public:
  // EMparticle(): particle(){};
  EMparticle(Double_t m0_ = 1,
             Double_t q_  = 1,
             TVector3 x0_ = TVector3(),
             TVector3 v0_ = TVector3(),
             Bool_t   Invincible_ = false,
             Bool_t   Record_     = false):
             particle(m0_, x0_, v0_, Invincible_, Record_){
               q = q_;
             }
  ~EMparticle(){}
  Double_t GetQ(){return q;}


private:
  // Double_t q;
};

#endif
