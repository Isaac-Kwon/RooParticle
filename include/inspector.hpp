#ifndef H_INSPECTOR
#define H_INSPECTOR 1;

#include "TROOT.h"
#include "particle.hpp"
#include "event.hpp"

class event;

class inspector{
  public:
  inspector(){}
  ~inspector(){}
  
  virtual Bool_t  Inspect(){return kFALSE;}
  virtual Float_t Evaluate(){return 0.;}
          void    SetValue(Float_t val_){val=val_;}
          void    SetMethod(TString methodname){ftn=methodname;}
          TString GetMethod(){return ftn;}
          Float_t GetCriterion(){return val;}
  virtual TString Print(Bool_t element=kFALSE, Bool_t mechanics=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);
  virtual Bool_t  IsInspectorP(){return kFALSE;}
  virtual Bool_t  IsInspectorE(){return kFALSE;}
  
  protected:
  TString ftn="";
  Float_t val=-1;

  private:
  // inspectorT(particle *p1_, particle *p2_, TString method_, Float_t val_=-1){};
  // inspectorT(event* event_, TString method_, Float_t val_=-1){};

};

class inspectorP: public inspector{ //Calculator, Static
  public:
  inspectorP(particle *p1_, particle *p2_, TString method_, Float_t val_=-1);
  
  Bool_t  Inspect();
  Float_t Evaluate(); //Get Evaluated value from inspecting method

  void SetParticles(particle * p1_, particle * p2_){p1 = p1_; p2 = p2_;}
  void SetParticle1(particle * p1_){p1 = p1_;}
  void SetParticle2(particle * p2_){p2 = p2_;}
  
  particle* GetParticle1(){return p1;}
  particle* GetParticle2(){return p2;}

  TString Print(Bool_t particles=kFALSE, Bool_t mechanics=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);

  // Bool_t SetInitial(Bool_t index=0, TString mtd="", Double_t var=0.1);
  Bool_t SetInitial_D(Bool_t index=0, TString mtd="", Double_t var=0.1);

  Bool_t IsInspectorP(){return kTRUE;}
  // "COS" : Cosine value 
  // "DEG" : Angle (Degree)
  // "RAD" : Angle (Radian)
  // "CNT" : DerivingCount

  protected:
  Float_t CalculateCosine();
  Float_t CalculateAngle(Bool_t rad=kFALSE);
  Long_t  CountDeriving();
  Float_t CaculateSPD();

  private:
  particle *p1, *p2;
};


class inspectorE: public inspector{
  public:
  inspectorE(event* event_, TString method_, Float_t val_=-1);

  void     SetEvent(event* event_){evt=event_;}
  event *  GetEvent(){return evt;}
  Bool_t   Inspect();
  Float_t Evaluate();

  Bool_t   IsInspectorE(){return kTRUE;}

  protected:

  private:
  event *evt;
  // TString ftn=TString("");
  // Float_t val=-1;

};

#endif