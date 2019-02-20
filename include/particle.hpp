
#ifndef H_PARTICLE
#define H_PARTICLE 1;

#include "TROOT.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "TVectorT.h"

class particle{
public:
  // particle();
  particle(Double_t m0_ = 1,
           //Double_t q_  = 1,
           TVectorD x0_ = TVectorD(3),
           TVectorD v0_ = TVectorD(3),
           Bool_t   Invincible_ = false);
  //Force Methods (Holding and Release)
  particle(const particle & p); //Copy Constructor
  ~particle(){delete ppath;} //Destructor (NO Pointer Variable in class)

  void holdForce(const TVectorD f);
  void releaseForce(const Double_t dt);
  void resetHold();

  //Applying Methods
  void applyMechanic(const TVectorD x_, const TVectorD v_){x = x_; v = v_;}
  void applyX(const TVectorD x_){x=x_;}
  void applyV(const TVectorD v_){v=v_;}
  void applyForce(const TVectorD f, const Double_t dt=1.);
  void applyDX(const TVectorD dx){x += dx;}
  void applyDV(const TVectorD dv){v += dv;}
  void applyDelta(const TVectorD dx, const TVectorD dv){x += dx; v += dv;}

  void SetName(const TString name_){name=name_;}

  TVectorD GetX(){return TVectorD(x);}
  TVectorD GetV(){return TVectorD(v);}
  TVectorD GetP(){return TVectorD(v * GetM());}
  Double_t GetM0(){return m0;}
  Double_t GetM(Bool_t restmass=false);

  Bool_t IsInvincible(){return Invincible;}

protected:
private:
  TString name;
  TVectorD x          = TVectorD(3);
  TVectorD v          = TVectorD(3);
  Double_t m0         = 1;
  //Double_t q          = 1;
  Bool_t   Invincible = false;

  Bool_t   hold       = false;
  TVectorD hF         = TVectorD(3);
};

class EMparticle : public particle{
public:
  // EMparticle(): particle(){};
  EMparticle(Double_t m0_ = 1,
             Double_t q_  = 1,
             TVectorD x0_ = TVectorD(3),
             TVectorD v0_ = TVectorD(3),
             Bool_t   Invincible_ = false):
             particle(m0_, x0_, v0_, Invincible_){
               q = q_;
             }
  Double_t GetQ(){return q;}

private:
  Double_t q;
};

#endif
