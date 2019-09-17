#include "iostream"

#include "TROOT.h"
#include "TString.h"
#include "TVectorT.h"
#include "TMath.h"

#include "inspector.hpp"

Float_t XVCosine(particle* p1, particle* p2){
    TVector3 dX = (p1->GetX())-(p2->GetX());
    TVector3 dV = (p1->GetV())-(p2->GetV());
    Double_t cc = (dX*dV)/(dX.Mag() * dV.Mag());
    return cc;
}

// Float_t XVCosine(TVector3 dX, TVector3 dV){
//   return (dX*dV)/(TMath::Sqrt(dX.Norm2Sqr() * dV.Norm2Sqr()));
// }

inspector::inspector(){
}

inspector::inspector(particle *p1_, particle *p2_, TString method_, Float_t val_){
    p1 = p1_;
    p2 = p2_;
    ftn = method_;
    val = val_;
}

inspector::~inspector(){
}

Float_t inspector::Evaluate(){
  if(ftn=="COS"){
    return CalculateCosine();
  }else if(ftn=="DEG"){
    return CalculateAngle(kFALSE);
  }else if(ftn=="RAD"){
    return CalculateAngle(kTRUE);
  }else{
    std::cout<<"No Function Name Defined / "<< ftn << std::endl;
    return kFALSE;
  }
  std::cout<<"INSPECTING METHOD IS NOT DEFINED"<<std::endl;
  return kFALSE;
}

Bool_t inspector::Inspect(){
    return Evaluate()<val;
}

Float_t inspector::CalculateCosine(){
    TVector3 dX = (p1->GetX())-(p2->GetX());
    TVector3 dV = (p1->GetV())-(p2->GetV());
    Double_t cc = (dX*dV)/(dX.Mag() * dV.Mag());
    return cc;
}

Float_t inspector::CalculateAngle(Bool_t rad){
    Double_t cc = CalculateCosine();
    if(rad){
      return TMath::ACos(cc);
    }else{
      return TMath::ACos(cc) * TMath::RadToDeg();
    }
}

Bool_t inspector::SetInitial_D(Bool_t index, TString mtd,Double_t var){ //set primary condition from present condition
  // if index=0(false) -> move p1, index=1(true) -> move p2
  
  // Initial position solver
  //    with dot product
  // Go around signal (x, v)
  //   angle between x and v -> less than epsilon

  // Approaching signal 
  //   angle between x and v -> larger than 1 (or pi) - eps
  cerr<<"inspector::SetInitial_D function needs to be fixed. do not use without any instruction "<<endl;

  if(mtd==""){
    mtd = ftn;
  }


  //DEFAULT, BASED ON COS_VALUE
  Double_t objval;
  if(mtd=="COS"){
    objval = -1+var;
  }else if(mtd=="DEG"){
    objval = TMath::Cos((180-var)*TMath::DegToRad());
  }else if(mtd=="RAD"){
    objval = TMath::Cos((180-var));
  }else{
    std::cerr<<"No Function Name Defined / "<< ftn << std::endl;
    return kFALSE;
  }

  particle *p0_, *p1_;
  if(index){
    p0_ = p1;
    p1_ = p2;
  }else{
    p0_ = p2;
    p1_ = p1;
  }

  Double_t cos0, cos1;
  cos0 = XVCosine(p0_, p1_);
  cos1 = objval;
  Double_t A1 = (p0_->GetX().Mag());///(p0_->GetV().Norm2Sqr()));
  Double_t A2 = cos0-TMath::Sqrt(1-TMath::Power(cos0,2))*(cos1/TMath::Sqrt(1-TMath::Power(cos1,2)));  
  Double_t A = A1*A2;

  p0_->x -= A*p0_->v;


  return kTRUE;

  // A = (|x0|/|v0|)*(cos(th0)-cos(th1)*sin(th0)*sqrt(1-(cos(th1))^2))
  // cos th1 :: settedvalue
}

Bool_t inspector::SetInitial(Bool_t index, TString mtd, Double_t var){
  return kTRUE;
}



TString inspector::Print(Bool_t particles, Bool_t mechanics, Bool_t mute, Bool_t pprint){
  TString result = TString("");
  result += TString("[ Inspector Particle: ");
  if(particles){
    result += GetParticle1()->Print(kTRUE, kTRUE);
    result += GetParticle2()->Print(kTRUE, kTRUE);
  }

  result += ", MTD: ";
  result += ftn;

  result += ", CRT: ";
  result += GetCriterion();

  result += ", VAL: ";
  result += Evaluate();

  result += " ]";
  
  if(!mute){
    cout<<result<<endl;
  }
  return result;
    
}


