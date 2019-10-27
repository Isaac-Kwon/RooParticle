#include "TROOT.h"
#include "TVector3.h"

#include "path.hpp"
#include "list"
#include "iostream"

#include <math.h>

using namespace std;

path::path(){
}

path::path(path& path_){
  std::vector<TVector3>::iterator i1;
  std::vector<Double_t>::iterator i2;

  for(i1=path_.pos.begin();i1!=path_.pos.end();++i1){
    pos.push_back(*i1);
  }

  for(i1=path_.vel.begin();i1!=path_.vel.end();++i1){
    vel.push_back(*i1);
  }

    for(i2=path_.t.begin();i2!=path_.t.end();++i2){
    t.push_back(*i2);
  }
}

path::~path(){
  pos.clear();
  vel.clear();
  t.clear();

  // pos.shrink_to_fit();
  // vel.shrink_to_fit();
  // t.shrink_to_fit();
}

void path::AddRecord(TVector3 x_, TVector3 v_, Double_t time_=-1.){
  pos.push_back(x_);
  vel.push_back(v_);
  t.push_back(time_);
}

TVector3 path::GetRecordX(Int_t i){
  // std::cout<<"CALL RECORDX ("<<i<<")"<<std::endl;
  if (pos.size()==0 || pos.size()<=i){
    std::cerr<<"Path::GetRecordX : NO ANY VALUES, RETURN NULL VECTOR"<<std::endl;
    return TVector3();
  }
  return pos[i];
}

TVector3 path::GetRecordV(Int_t i){
  if (vel.size()==0 || vel.size()<=i){
    std::cerr<<"Path::GetRecordV : NO ANY VALUES, RETURN NULL VECTOR"<<std::endl;
    return TVector3();
  }
  return vel[i];
}

Double_t path::GetRecordT(Int_t i){
  if (t.size()==0 || t.size()<=i){
    std::cerr<<"Path::GetRecordT : NO ANY VALUES, RETURN NULL VALUE"<<std::endl;
    return Double_t(0.);
  }
  return t[i];
}

void path::GetRecord(Int_t i,
                     TVector3 & x_, TVector3 & v_,
                     Double_t & t_){
  x_ = GetRecordX(i);
  v_ = GetRecordV(i);
  t_ = GetRecordT(i);
  return;
}

TVector3 path::GetInitX(){
  if (pos.size()==0){
    std::cerr<<"Path::GetInitX : RETURN NULL VECTOR"<<std::endl;
    return TVector3();
  }
  return pos.front();
}

TVector3 path::GetInitV(){
  if (vel.size()==0){
    std::cerr<<"Path::GetInitV : NO ANY VALUES, RETURN NULL VECTOR"<<std::endl;
    return TVector3();
  }
  return vel.front();
}

Double_t path::GetInitT(){
  if (t.size()==0){
    std::cerr<<"Path::GetInitT : NO ANY VALUES, RETURN ZERO VALUE"<<std::endl;
    return Double_t(0.);
  }
  return t.front();
}

void path::GetInit(TVector3 & x_, TVector3 & v_,
                   Double_t & t_){
  x_ = GetInitX();
  v_ = GetInitV();
  t_ = GetInitT();
  return;
}

TVector3 path::GetLastX(){
  if (pos.size()==0){
    std::cerr<<"Path::GetLastX : NO ANY VALUES, RETURN NULL VECTOR"<<std::endl;
    return TVector3();
  }
  return pos.back();
}

TVector3 path::GetLastV(){
  if (vel.size()==0){
    std::cerr<<"Path::GetLastV : NO ANY VALUES, RETURN NULL VECTOR"<<std::endl;
    return TVector3();
  }
  return vel.back();
}

Double_t path::GetLastT(){
  if (t.size()==0){
    std::cerr<<"Path::GetLastT : NO ANY VALUES, RETURN ZERO VALUE"<<std::endl;
    return Double_t(0.);
  }
  return t.back();
}

void path::GetLast(TVector3 & x_, TVector3 & v_,
                   Double_t & t_){
  x_ = GetLastX();
  v_ = GetLastV();
  t_ = GetLastT();
  return;
}

void path::GetDCA(const TVector3 pos_, Double_t & DCA, Double_t & outtime, Bool_t verbose){
  Int_t i;
  Double_t resultDCA=std::numeric_limits<double>::max();
  Double_t tempDCA=-2.;
  Int_t resultTI=-1;
  TVector3 tempvector;
  // for(std::vector<TVector3>::iterator it = pos.begin();
  //     it != pos.end(); ++it){
  for(i=0; i<pos.size(); i++){
    tempvector = pos.at(i)-GetRecordX(i);
    tempDCA = tempvector.Mag();
    if(tempDCA<resultDCA){
      resultDCA = tempDCA;
      resultTI = i;
    }
  }
  if(resultTI==-1 || verbose){
    std::cerr<<"path::GetDCA : Algorithm Error - Result Index 0"<<std::endl;
  }
  DCA     = resultDCA;
  outtime = t[resultTI];
}

void path::GetDCA(path* timepath, Double_t & DCA, Double_t & outtime, Bool_t verbose){
  if(!timeChk(timepath)){
    return;
  }

  Double_t resultDCA=std::numeric_limits<double>::max();
  Double_t tempDCA=-2.;
  Int_t resultTI=-1;
  TVector3 tempvector;

  Int_t i;
  for(i=0; i<pos.size(); ++i){
    tempvector = (timepath->GetRecordX(i))-GetRecordX(i);
    tempDCA = tempvector.Mag();
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
    std::cerr<<"path::lengthChk : LENGTH CHK ERROR (NOT SAME)"<<std::endl;
    return kFALSE;
  }
  else{
    return kTRUE;
  }
}

Bool_t path::timeChk(path* target){
  if(!lengthChk(target)){ 
    return kFALSE;
  }
  Int_t i;
  for(i=0; i<target->GetMaxNumber(); i++){
    if(target->GetRecordT(i)!=GetRecordT(i)){
      std::cerr<<"path::timeChk : TIME CHK ERROR: "<< i << std::endl;
      return false;
    }
  }
  return true;
}
