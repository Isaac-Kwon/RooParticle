#ifndef H_PATH
#define H_PATH 1;

#include "TROOT.h"
#include "RPConfig.hpp"

class path{
public:
  ~path(){};
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
protected:
private:
  Int_t writingCursor = 0;

  Double_t t[MAXSTEPS];

  Double_t x[MAXSTEPS];
  Double_t y[MAXSTEPS];
  Double_t z[MAXSTEPS];

  Double_t vx[MAXSTEPS];
  Double_t vy[MAXSTEPS];
  Double_t vz[MAXSTEPS];
};

#endif
