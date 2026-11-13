#include "iostream"

#include "TROOT.h"
#include "TString.h"
#include "TVectorT.h"
#include "TMath.h"

#include "inspector.hpp"
#include "asymptoteN.hpp"

#include <limits>

Double_t XVCosine(particle* p1, particle* p2){
    TVector3 dX = (p1->GetX())-(p2->GetX());
    TVector3 dV = (p1->GetV())-(p2->GetV());
    Double_t cc = (dX*dV)/(dX.Mag() * dV.Mag());
    return cc;
}

// Double_t XVCosine(TVector3 dX, TVector3 dV){
//   return (dX*dV)/(TMath::Sqrt(dX.Norm2Sqr() * dV.Norm2Sqr()));
// }

// The name is parsed once (SetMethod) ; names and enumerators are listed here only.
static const struct {const char * name; inspector::methodID id;} gMethodTable[] = {
  {"COS",  inspector::kCOS},  {"DEG",  inspector::kDEG},  {"RAD",  inspector::kRAD},
  {"CNT",  inspector::kCNT},  {"SPDG", inspector::kSPDG}, {"SPDL", inspector::kSPDL},
  {"DSTG", inspector::kDSTG}, {"DSTL", inspector::kDSTL}, {"VRG",  inspector::kVRG},
  {"VRL",  inspector::kVRL},
  {"KEG",  inspector::kKEG},  {"KEL",  inspector::kKEL},  {"PEG",  inspector::kPEG},
  {"PEL",  inspector::kPEL},  {"PKG",  inspector::kPKG},  {"PKL",  inspector::kPKL},
  {"EDG",  inspector::kEDG},  {"EDL",  inspector::kEDL},  {"TMG",  inspector::kTMG},
  {"TML",  inspector::kTML},  {"DXG",  inspector::kDXG},  {"DXL",  inspector::kDXL},
  {"DNG",  inspector::kDNG},  {"DNL",  inspector::kDNL},  {"DSG",  inspector::kDSG},
  {"DSL",  inspector::kDSL},  {"ACG",  inspector::kACG},  {"ACL",  inspector::kACL},
  {"VOUT", inspector::kVOUT}, {"VOUTA",inspector::kVOUTA},
  {"ASYD", inspector::kASYD}, {"ASYV", inspector::kASYV}
};

inspector::methodID inspector::ParseMethod(TString name){
  for(const auto & e : gMethodTable){
    if(name==e.name) return e.id;
  }
  return kInvalid;
}

TString inspector::MethodName(methodID id){
  for(const auto & e : gMethodTable){
    if(e.id==id) return TString(e.name);
  }
  return TString("");
}

inspectorP::inspectorP(particle *p1_, particle *p2_, TString method_, Double_t val_): inspector(){
  SetParticles(p1_,p2_);
  SetMethod(method_);
  SetValue(val_);
}

inspectorP::inspectorP(particle *p1_, particle *p2_, methodID method_, Double_t val_): inspector(){
  SetParticles(p1_,p2_);
  SetMethod(method_);
  SetValue(val_);
}

Double_t inspectorP::Evaluate(){
  switch(mtd){
    case kCOS:  return CalculateCosine();
    case kDEG:  return CalculateAngle(kFALSE);
    case kRAD:  return CalculateAngle(kTRUE);
    case kCNT:  return CountDeriving();
    case kSPDG: case kSPDL: return CaculateSPD();
    case kDSTG: case kDSTL: return CalculateDistance();
    case kVRG:  case kVRL:  return CalculateRadialVelocity();
    default:
      std::cout<<"No Function Name Defined / "<< ftn << std::endl;
      return kFALSE;
  }
}

Bool_t inspectorP::IsValid(){
  return mtd>=kCOS && mtd<=kVRL;
}

Bool_t inspectorP::Inspect(){
    if(!IsValid()){
      return kFALSE;
    }
    // Undefined angle / cosine / radial velocity : never satisfied.
    // (not left to NaN comparisons ; atan2(0,0) would give 0 degree and satisfy "<")
    if((mtd==kCOS || mtd==kDEG || mtd==kRAD) && IsDegenerate(kTRUE)){
      return kFALSE;
    }
    if((mtd==kVRG || mtd==kVRL) && IsDegenerate(kFALSE)){
      return kFALSE;
    }
    if(mtd==kCNT || mtd==kSPDG || mtd==kDSTG || mtd==kVRG){
      return Evaluate()>val;
    }
    return Evaluate()<val;
}

Bool_t inspectorP::IsDegenerate(Bool_t needV){
  // exact zero : the direction does not exist
  if(((p1->GetX())-(p2->GetX())).Mag2()==0.) return kTRUE;
  if(needV && ((p1->GetV())-(p2->GetV())).Mag2()==0.) return kTRUE;
  return kFALSE;
}

Double_t inspectorP::CalculateDistance(){
  return ((p1->GetX())-(p2->GetX())).Mag();
}

Double_t inspectorP::CalculateRadialVelocity(){
  // positive : the two particles are receding from each other
  TVector3 dX = (p1->GetX())-(p2->GetX());
  TVector3 dV = (p1->GetV())-(p2->GetV());
  if(dX.Mag2()==0.) return TMath::QuietNaN();
  return (dX*dV)/dX.Mag();
}

Double_t inspectorP::CalculateCosine(){
    if(IsDegenerate(kTRUE)) return TMath::QuietNaN();
    TVector3 dX = (p1->GetX())-(p2->GetX());
    TVector3 dV = (p1->GetV())-(p2->GetV());
    Double_t cc = (dX*dV)/(dX.Mag() * dV.Mag());
    // rounding can push it slightly out of [-1, 1]
    return TMath::Max(-1., TMath::Min(1., cc));
}

Double_t inspectorP::CalculateAngle(Bool_t rad){
    if(IsDegenerate(kTRUE)) return TMath::QuietNaN();
    // atan2(|dX x dV|, dX.dV) keeps the precision near 0 and pi, where acos(cos) does not
    TVector3 dX = (p1->GetX())-(p2->GetX());
    TVector3 dV = (p1->GetV())-(p2->GetV());
    Double_t angle = TMath::ATan2(dX.Cross(dV).Mag(), dX.Dot(dV));
    return rad ? angle : angle * TMath::RadToDeg();
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

inspectorE::inspectorE(event* event_, methodID method_, Double_t val_): inspector(){
  SetEvent(event_);
  SetMethod(method_);
  SetValue(val_);
}

Bool_t inspectorE::IsValid(){
  return mtd>=kKEG && mtd<=kASYV;
}

Bool_t inspectorE::Inspect(){
  switch(mtd){
    case kPEG: return (evt->GetNetPE()>val);
    case kPEL: return (evt->GetNetPE()<val);
    case kKEG: case kPKG: case kEDG: case kTMG: case kDSG: case kACG:
      return (Evaluate()>val);
    case kKEL: case kPKL: case kEDL: case kTML: case kDSL: case kACL:
      return (Evaluate()<val);
    case kVOUT:  return !(evt->IsInsideVolume(kTRUE));   // not all inside : somebody is out
    case kVOUTA: return !(evt->IsInsideVolume(kFALSE));  // nobody inside : all are out
    case kASYD: case kASYV: {
      // The estimate only means "what is left" once every interacting pair
      // is receding : before that the collision is still ahead.
      Bool_t separating = kFALSE;
      Double_t rest = (mtd==kASYD) ? asymptoteN::RemainingDeflection(evt, &separating)
                                   : asymptoteN::RemainingChange(evt, &separating);
      return separating && (rest<val);
    }
    case kDXG: return (evt->GetDeriveN(kFALSE)>val);
    case kDXL: return (evt->GetDeriveN(kTRUE)<val);
    case kDNG: return (evt->GetDeriveN(kFALSE)>val);
    case kDNL: return (evt->GetDeriveN(kTRUE)<val);
    default:
      std::cerr<<"inspectorE::Inspect : inspecting method isn't defined well." << std::endl;
      return kFALSE;
  }
}

Double_t inspectorE::Evaluate(){
  switch(mtd){
    // Kinetic energy : the relativistic one when the relativistic mode is
    // really in effect (Taylor ignores the setting), as for "EDG"/"EDL".
    case kKEG: case kKEL:
      return (evt->IsRelativistic() && evt->GetIntegrator()!=event::kTaylor) ? evt->GetNetKERel() : evt->GetNetKE();
    case kPEG: case kPEL: return evt->GetNetPE();
    case kDXG: case kDXL: return evt->GetDeriveN(kFALSE);
    case kDNG: case kDNL: return evt->GetDeriveN(kTRUE);
    case kPKG: case kPKL: {
      Double_t ke = (evt->IsRelativistic() && evt->GetIntegrator()!=event::kTaylor) ? evt->GetNetKERel() : evt->GetNetKE();
      if(ke==0.){
        return std::numeric_limits<Double_t>::max();
      }
      return TMath::Abs(evt->GetNetPE()/ke);
    }
    case kEDG: case kEDL: {
      // In relativistic mode only the relativistic energy is conserved.
      // (Taylor ignores the relativistic setting.)
      Bool_t rel = evt->IsRelativistic() && evt->GetIntegrator()!=event::kTaylor;
      Double_t e = rel ? evt->GetNetERel() : evt->GetNetE();
      if(!hasRef){
        ref = e;
        hasRef = kTRUE;
      }
      if(ref==0.){
        return TMath::Abs(e);
      }
      return TMath::Abs((e-ref)/ref);
    }
    case kTMG: case kTML: return evt->GetTime();
    case kDSG: case kDSL: return evt->GetMinDistance();
    case kACG: case kACL: return evt->GetMaxAcceleration();
    case kVOUT:  return !(evt->IsInsideVolume(kTRUE));
    case kVOUTA: return !(evt->IsInsideVolume(kFALSE));
    case kASYD: return asymptoteN::RemainingDeflection(evt);
    case kASYV: return asymptoteN::RemainingChange(evt);
    default:
      std::cerr<<"inspectorE::Inspect : inspecting method isn't defined well." << std::endl;
      return kFALSE;
  }
}

TString inspectorE::Print(Bool_t element, Bool_t mechanics, Bool_t mute, Bool_t pprint){
  TString result = TString("");
  result += TString("[ Inspector Event: ");

  result += "MTD: ";
  result += ftn;

  result += ", CRT: ";
  result += GetCriterion();

  result += ", VAL: ";
  result += Evaluate();

  result += ", RST: ";
  result += Inspect();

  if(IsAbort()){
    result += ", ABORT";
  }

  result += " ]";

  if(!mute){
    cout<<result<<endl;
  }
  return result;
}


