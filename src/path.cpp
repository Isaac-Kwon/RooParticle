#include "TROOT.h"
#include "TVectorD.h"

#include "path.hpp"
#include "iostream"

using namespace std;

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
