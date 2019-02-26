#ifndef H_PATH
#define H_PATH 1;

#include "TROOT.h"
#include "RPConfig.hpp"

class path{
public:
  ~path();
  Int_t AddRecord(TVectorD x_, TVectorD v_, Double_t time=-1.);
  TVectorD GetRecordX(Int_t i);
  TVectorD GetRecordV(Int_t i);
  TVectorD GetInitX();
  TVectorD GetInitV();
  TVectorD GetLastX();
  TVectorD GetLastV();
  Int_t    GetMaxNumber(){return writingCursor;};
  Double_t * GetAllX(){return x;}
  Double_t * GetAllY(){return y;}
  Double_t * GetAllZ(){return z;}
  Double_t * GetAllVX(){return vx;}
  Double_t * GetAllVY(){return vy;}
  Double_t * GetAllVZ(){return vz;}
  void     GetDCA(const TVectorD pos, Double_t & DCA, Double_t & outtime, Bool_t verbose=false);
  void     GetDCA(path* timepath, Double_t & DCA, Double_t & outtime, Bool_t verbose=false);
  // Double_t GetDCA(const TVectorD pos);
  // Double_t GetDCA(path* timepath);
protected:
private:
  Int_t writingCursor = 0;

  Double_t * t = new Double_t[MAXSTEPS];

  Double_t * x = new Double_t[MAXSTEPS];
  Double_t * y = new Double_t[MAXSTEPS];
  Double_t * z = new Double_t[MAXSTEPS];

  Double_t * vx = new Double_t[MAXSTEPS];
  Double_t * vy = new Double_t[MAXSTEPS];
  Double_t * vz = new Double_t[MAXSTEPS];

  Bool_t lengthChk(path* target);
  Bool_t timeChk(path* target);
};

#endif
