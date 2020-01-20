#include "iostream"

#include "TROOT.h"
#include "TString.h"
#include "TVectorT.h"
#include "TMath.h"

#include "inspector.hpp"

Double_t XVCosine(particle* p1, particle* p2){
    TVector3 dX = (p1->GetX())-(p2->GetX());
    TVector3 dV = (p1->GetV())-(p2->GetV());
    Double_t cc = (dX*dV)/(dX.Mag() * dV.Mag());
    return cc;
}

// Double_t XVCosine(TVector3 dX, TVector3 dV){
//   return (dX*dV)/(TMath::Sqrt(dX.Norm2Sqr() * dV.Norm2Sqr()));
// }

inspectorP::inspectorP(particle *p1_, particle *p2_, TString method_, Double_t val_): inspector(){
  SetParticles(p1_,p2_);
  SetMethod(method_);
  SetValue(val_);
}

Double_t inspectorP::Evaluate(){
  if(ftn=="COS"){
    return CalculateCosine();
  }else if(ftn=="DEG"){
    return CalculateAngle(kFALSE);
  }else if(ftn=="RAD"){
    return CalculateAngle(kTRUE);
  }else if(ftn=="CNT"){
    return CountDeriving();
  }else if(ftn=="SPDG" || ftn=="SPDL"){
    return CaculateSPD();
  }else{
    std::cout<<"No Function Name Defined / "<< ftn << std::endl;
    return kFALSE;
  }
  std::cout<<"INSPECTING METHOD IS NOT DEFINED"<<std::endl;
  return kFALSE;
}

Bool_t inspectorP::Inspect(){
    if(ftn=="CNT" || ftn=="SPDG"){
      return Evaluate()>val;
    }
    return Evaluate()<val;
}

Double_t inspectorP::CalculateCosine(){
    TVector3 dX = (p1->GetX())-(p2->GetX());
    TVector3 dV = (p1->GetV())-(p2->GetV());
    Double_t cc = (dX*dV)/(dX.Mag() * dV.Mag());
    return cc;
}

Double_t inspectorP::CalculateAngle(Bool_t rad){
    Double_t cc = CalculateCosine();
    if(rad){
      return TMath::ACos(cc);
    }else{
      return TMath::ACos(cc) * TMath::RadToDeg();
    }
}

Double_t inspectorP::CaculateSPD(){
  return ((p1->GetV())-(p2->GetV())).Mag();
}

Long_t inspectorP::CountDeriving(){
  return min(p1->GetECNT(), p2->GetECNT());

}

Bool_t inspectorP::SetInitial_D(Bool_t index, TString mtd,Double_t var){ //set primary condition from present condition
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

  p0_->x -= A*(p0_->GetV());


  return kTRUE;

  // A = (|x0|/|v0|)*(cos(th0)-cos(th1)*sin(th0)*sqrt(1-(cos(th1))^2))
  // cos th1 :: settedvalue
}

// Bool_t inspectorP::SetInitial(Bool_t index, TString mtd, Double_t var){
//   return kTRUE;
// }

TString inspector::Print(Bool_t element, Bool_t mechanics, Bool_t mute, Bool_t pprint){
  return TString("");
}



TString inspectorP::Print(Bool_t particles, Bool_t mechanics, Bool_t mute, Bool_t pprint){
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

  result += ", RST: ";
  result += Inspect();

  result += " ]";
  
  if(!mute){
    cout<<result<<endl;
  }
  return result;
    
}


//========================================//

inspectorE::inspectorE(event* event_, TString method_, Double_t val_): inspector(){
  SetEvent(event_);
  SetMethod(method_);
  SetValue(val_);
}

Bool_t inspectorE::Inspect(){
  if(ftn=="KEG"){
    return (evt->GetNetKE()>val);
  }else if(ftn=="KEL"){
    return (evt->GetNetKE()<val);
  }else if(ftn=="PEG"){
    return (evt->GetNetPE()>val);
  }else if(ftn=="PEL"){
    return (evt->GetNetPE()<val);
  }else if(ftn=="DXG"){
    return (evt->GetDeriveN(kFALSE)>val);
  }else if(ftn=="DXL"){
    return (evt->GetDeriveN(kTRUE)<val);
  }else if(ftn=="DNG"){
    return (evt->GetDeriveN(kFALSE)>val);  
  }else if(ftn=="DNL"){
    return (evt->GetDeriveN(kTRUE)<val);  
  }else{
    std::cerr<<"inspectorE::Inspect : inspecting method isn't defined well." << std::endl;
    return kFALSE;
  }
  return kFALSE;
}

Double_t inspectorE::Evaluate(){
  if(ftn=="KEG" || ftn=="KEL"){
    return evt->GetNetKE();
  }else if(ftn=="PEG" || ftn=="PEL"){
    return evt->GetNetPE();
  }else if(ftn=="DXG" || ftn=="DXL"){
    return evt->GetDeriveN(kFALSE);
  }else if(ftn=="DNG" || ftn=="DNL"){
    return evt->GetDeriveN(kTRUE);  
  }else{
    std::cerr<<"inspectorE::Inspect : inspecting method isn't defined well." << std::endl;
    return kFALSE;
  }
  return kFALSE;
}


