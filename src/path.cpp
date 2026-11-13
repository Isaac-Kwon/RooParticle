#include "TROOT.h"
#include "TVector3.h"

#include "path.hpp"
#include "list"
#include "iostream"

#include <math.h>

using namespace std;

path::path(){
}

path::path(const path& path_){
  std::vector<TVector3>::const_iterator i1;
  std::vector<Double_t>::const_iterator i2;

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

// Sentinel written to DCA and outtime when no distance can be computed
// (empty path, length/time mismatch, null pointer).
static const Double_t kPathDCAInvalid = -1.;

// Refine a discrete minimum of the squared distance d2[] at index k by a
// parabola through samples k-1, k, k+1 (t need not be uniform). The recorded
// steps are coarse, so the true closest approach usually lies between two
// of them. Keeps the discrete values at the path ends, for non-increasing
// time, or when the parabola has no interior minimum (degenerate fit).
static void RefineMinimum(const std::vector<Double_t> & d2,
                          const std::vector<Double_t> & t,
                          Int_t k, Double_t & DCA, Double_t & outtime){
  DCA     = sqrt(d2[k]);
  outtime = t[k];
  if(k<=0 || k>=(Int_t)d2.size()-1) return;

  Double_t s0 = t[k-1]-t[k];
  Double_t s2 = t[k+1]-t[k];
  if(!(s0<0. && s2>0.)) return;

  Double_t y1 = d2[k];
  Double_t q0 = (d2[k-1]-y1)/s0;
  Double_t q2 = (d2[k+1]-y1)/s2;
  Double_t a  = (q0-q2)/(s0-s2);
  if(!(a>0.)) return;
  Double_t b  = q2-a*s2;
  Double_t sm = -b/(2.*a);
  if(sm<s0 || sm>s2) return;

  Double_t dmin2 = y1-b*b/(4.*a);
  if(dmin2<0.) dmin2 = 0.; // rounding only; distance squared cannot be negative
  DCA     = sqrt(dmin2);
  outtime = t[k]+sm;
}

void path::GetDCA(const TVector3 pos_, Double_t & DCA, Double_t & outtime, Bool_t verbose, Bool_t interpolate){
  DCA     = kPathDCAInvalid;
  outtime = kPathDCAInvalid;

  if(pos.empty() || pos.size()!=t.size()){
    std::cerr<<"path::GetDCA : EMPTY OR INCONSISTENT PATH, RETURN DCA=-1, TIME=-1"<<std::endl;
    return;
  }

  std::vector<Double_t> d2(pos.size());
  Int_t resultTI=0;
  for(Int_t i=0; i<(Int_t)pos.size(); i++){
    // distance to the fixed point (not to the path itself)
    d2[i] = (pos[i]-pos_).Mag2();
    if(d2[i]<d2[resultTI]) resultTI = i;
  }

  if(interpolate){
    RefineMinimum(d2, t, resultTI, DCA, outtime);
  }
  else{
    DCA     = sqrt(d2[resultTI]);
    outtime = t[resultTI];
  }
  if(verbose){
    std::cerr<<"path::GetDCA : index "<<resultTI<<" DCA "<<DCA<<" time "<<outtime<<std::endl;
  }
}

void path::GetDCA(path* timepath, Double_t & DCA, Double_t & outtime, Bool_t verbose, Bool_t interpolate){
  DCA     = kPathDCAInvalid;
  outtime = kPathDCAInvalid;

  if(timepath==nullptr){
    std::cerr<<"path::GetDCA : NULL PATH, RETURN DCA=-1, TIME=-1"<<std::endl;
    return;
  }
  if(!timeChk(timepath)){
    std::cerr<<"path::GetDCA : PATHS DO NOT SHARE TIME GRID, RETURN DCA=-1, TIME=-1"<<std::endl;
    return;
  }
  if(pos.empty() || pos.size()!=t.size() || timepath->pos.size()!=pos.size()){
    std::cerr<<"path::GetDCA : EMPTY OR INCONSISTENT PATH, RETURN DCA=-1, TIME=-1"<<std::endl;
    return;
  }

  std::vector<Double_t> d2(pos.size());
  Int_t resultTI=0;
  for(Int_t i=0; i<(Int_t)pos.size(); ++i){
    d2[i] = (timepath->pos[i]-pos[i]).Mag2();
    if(d2[i]<d2[resultTI]) resultTI = i;
  }

  if(interpolate){
    RefineMinimum(d2, t, resultTI, DCA, outtime);
  }
  else{
    DCA     = sqrt(d2[resultTI]);
    outtime = t[resultTI];
  }
  if(verbose){
    std::cerr<<"path::GetDCA : index "<<resultTI<<" DCA "<<DCA<<" time "<<outtime<<std::endl;
  }
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
