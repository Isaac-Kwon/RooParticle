#ifndef H_INSPECTOR
#define H_INSPECTOR 1;

// #include "TROOT.h"
#include "particle.hpp"

using namespace std;

class inspector{ //Calculator, Static
  public:
  inspector();
  inspector(particle *p1_, particle *p2_, TString method_, Float_t val_=-1);
  virtual ~inspector();
  virtual Bool_t Inspect();

  void SetParticles(particle * p1_, particle * p2_){p1 = p1_; p2 = p2_;}
  void SetParticle1(particle * p1_){p1 = p1_;}
  void SetParticle2(particle * p2_){p2 = p2_;}
  void SetMethod(TString methodname){ftn=methodname;}

  particle* GetParticle1(){return p1;}
  particle* GetParticle2(){return p2;}

  TString GetMethod(){return ftn;}
  Float_t GetCriterion(){return val;}

  TString Print(Bool_t particles=kFALSE, Bool_t mechanics=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);

  Bool_t SetInitial(Bool_t index=0, TString mtd="", Double_t var=0.1);
  Bool_t SetInitial_D(Bool_t index=0, TString mtd="", Double_t var=0.1);
  
  // if index=0(false) -> move p1, index=1(true) -> move p2

  virtual Float_t Evaluate(); //Get Evaluated value from inspecting method

  // "COS" : Cosine value 
  // "DEG" : Angle (Degree)
  // "RAD" : Angle (Radian)
  // "CNT" : DerivingCount

  protected:
  Float_t CalculateCosine();
  Float_t CalculateAngle(Bool_t rad=kFALSE);
  Float_t CountDeriving();

  private:
  particle *p1, *p2;
  TString ftn="";
  Float_t val=-1;

  // TString dftn;
};

#endif