#ifndef H_PATH
#define H_PATH 1;

#include "TROOT.h"
#include "RPConfig.hpp"
#include "TVector3.h"
#include "iostream"

#include "list"

using namespace std;

class path{
public:
  path();
  path(path& path_);
  ~path();
  void AddRecord(TVector3 x_, TVector3 v_, Double_t time_);

  Int_t    GetRecordSize(){return pos.size();}

  TVector3 GetRecordX(Int_t i);
  TVector3 GetRecordV(Int_t i);
  Double_t GetRecordT(Int_t i);
  void     GetRecord(Int_t i, TVector3 & x_, TVector3 & v_, Double_t & t_);

  TVector3 GetInitX();
  TVector3 GetInitV();
  Double_t GetInitT();
  void     GetInit(TVector3 & x_, TVector3 & v_, Double_t & t_);

  TVector3 GetLastX();
  TVector3 GetLastV();
  Double_t GetLastT();
  void     GetLast(TVector3 & x_, TVector3 & v_, Double_t & t_);

  void     GetDCA(const TVector3 pos, Double_t & DCA, Double_t & outtime, Bool_t verbose=false);
  void     GetDCA(path* timepath, Double_t & DCA, Double_t & outtime, Bool_t verbose=false);
  Long64_t GetMaxNumber(){return t.size();}
  // Double_t GetDCA(const TVector3 pos);
  // Double_t GetDCA(path* timepath);
protected:
  Bool_t lengthChk(path* target);
  Bool_t timeChk(path* target);

private:
  // Int_t writingCursor = 0;
  std::vector<TVector3>  pos;
  std::vector<TVector3>  vel;
  std::vector<Double_t>  t;
};

#endif
