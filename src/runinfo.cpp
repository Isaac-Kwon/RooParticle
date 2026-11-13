#include "runinfo.hpp"

#include "TDirectory.h"
#include "TNamed.h"

#ifndef RP_GIT_SHA
#define RP_GIT_SHA "unknown"
#endif

UInt_t runinfo::SeedFromJob(Int_t jobIndex, Long64_t seed){
  if(seed >= 0) return (UInt_t) seed;
  return kDefaultSeed + (UInt_t) jobIndex;
}

TString runinfo::GitSHA(){
  return TString(RP_GIT_SHA);
}

void runinfo::Write(const char* key, const char* value){
  if(!gDirectory) return;
  TNamed info(key, value);
  // Write() keeps the key even though the stack object dies here.
  info.Write(key, TObject::kOverwrite);
}

void runinfo::Write(const char* key, Double_t value){
  Write(key, TString::Format("%.17g", value).Data());
}

void runinfo::Write(const char* key, Long64_t value){
  Write(key, TString::Format("%lld", value).Data());
}

void runinfo::WriteCommon(UInt_t seed, Double_t dt, const char* integrator,
                          Long64_t nevents, const char* inspectorCriteria){
  Write("seed", (Long64_t) seed);
  Write("dt", dt);
  Write("integrator", integrator);
  Write("nevents", nevents);
  Write("inspectors", inspectorCriteria);
  Write("git_sha", GitSHA().Data());
}
