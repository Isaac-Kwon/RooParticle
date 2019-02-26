#include "TROOT.h"
#include "TVectorD.h"

#include "path.hpp"
#include "iostream"

#include <math.h>

using namespace std;

path::~path(){
  delete[] t;

  delete[] x;
  delete[] y;
  delete[] z;

  delete[] vx;
  delete[] vy;
  delete[] vz;
}

Int_t path::AddRecord(TVectorD x_, TVectorD v_, Double_t time_){
  x[writingCursor] = x_[0];
  y[writingCursor] = x_[1];
  z[writingCursor] = x_[2];

  vx[writingCursor] = v_[0];
  vy[writingCursor] = v_[1];
  vz[writingCursor] = v_[2];

  t[writingCursor] = time_;

  writingCursor++;
  // std::cout<<"recorded"<<endl;
  return writingCursor-1;
}

TVectorD path::GetRecordX(Int_t i){
  if (writingCursor==0 || i>writingCursor){
    std::cout<<"NO ANY VALUES, RETURN NULL VECTOR"<<endl;
    return TVectorD(0);
  }
  Double_t result_[] = {x[i],y[i],z[i]};
  TVectorD result = TVectorD(3, result_);
  return result;
}

TVectorD path::GetRecordV(Int_t i){
  if (writingCursor==0 || i>writingCursor){
    std::cout<<"NO ANY VALUES, RETURN NULL VECTOR"<<endl;
    return TVectorD(0);
  }
  Double_t result_[] = {vx[i],vy[i],vz[i]};
  TVectorD result = TVectorD(3, result_);
  return result;
}

TVectorD path::GetInitX(){
  if (writingCursor==0){
    std::cout<<"NO ANY VALUES, RETURN NULL VECTOR"<<endl;
    return TVectorD(0);
  }
  return GetRecordX(0);
}

TVectorD path::GetInitV(){
  if (writingCursor==0){
    std::cout<<"NO ANY VALUES, RETURN NULL VECTOR"<<endl;
    return TVectorD(0);
  }
  return GetRecordV(0);
}

TVectorD path::GetLastX(){
  if (writingCursor==0){
    std::cout<<"NO ANY VALUES, RETURN NULL VECTOR"<<endl;
    return TVectorD(0);
  }
  return GetRecordX(writingCursor-1);
}

TVectorD path::GetLastV(){
  if (writingCursor==0){
    std::cout<<"NO ANY VALUES, RETURN NULL VECTOR"<<endl;
    return TVectorD(0);
  }
  return GetRecordV(writingCursor-1);
}

// TVectorD path::GetCloestPosition()

void path::GetDCA(const TVectorD pos, Double_t & DCA, Double_t & outtime, Bool_t verbose){
  Int_t i;
  Double_t resultDCA=std::numeric_limits<double>::max();
  Double_t tempDCA=-2.;
  Int_t resultTI=-1;
  TVectorD tempvector = TVectorD(3);
  for(i=0; i<writingCursor; i++){
    tempvector = pos-GetRecordX(i);
    tempDCA = sqrt(tempvector.Norm2Sqr());
    if(tempDCA<resultDCA){
      resultDCA = tempDCA;
      resultTI = i;
    }
  }
  if(resultTI==-1 || verbose){
    std::cout<<"Algorithm Error:: Result Index 0"<<endl;
  }
  DCA     = resultDCA;
  outtime = t[resultTI];
}

void path::GetDCA(path* timepath, Double_t & DCA, Double_t & outtime, Bool_t verbose){
  if(!timeChk(timepath)){
    return;
  }

  Int_t i;
  Double_t resultDCA=std::numeric_limits<double>::max();
  Double_t tempDCA=-2.;
  Int_t resultTI=-1;
  TVectorD tempvector = TVectorD(3);

  for(i=0; i<writingCursor; i++){
    tempvector = (timepath->GetRecordX(i))-GetRecordX(i);
    tempDCA = sqrt(tempvector.Norm2Sqr());
    if(tempDCA<resultDCA){
      resultDCA = tempDCA;
      resultTI = i;
    }
  }
  DCA     = resultDCA;
  outtime = t[resultTI];
}

Bool_t path::lengthChk(path* target){
  if(GetMaxNumber()!=target->GetMaxNumber()){
    std::cout<<"LENGTH CHK ERROR (NOT SAME)"<<endl;
    return false;
  }
  else{
    return true;
  }
}

Bool_t path::timeChk(path* target){
  if(!lengthChk(target)){
    return false;
  }

  Int_t i;
  for(i=0; i< target->GetMaxNumber();i++){
    if(target->t[i]!=t[i]){
      std::cout<<"TIME CHK ERROR: "<< i << endl;
      return false;
    }
  }
  return true;
}
