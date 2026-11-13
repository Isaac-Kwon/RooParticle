/// \file preEvent.hpp
/// \brief Unfinished, disabled feature: a "pre-event" that steps until two particles are aligned.
///
/// \warning Not implemented. The whole class preEvent (an event whose
/// DeriveInspect() steps until the angle between the relative position and
/// relative velocity of two particles falls below eps) is commented out here
/// and in src/preEvent.cpp; only an include of TROOT.h remains active.
/// Nothing in the code base uses it.
/// \todo Finish and port to the current event interface, or remove.

#ifndef H_PREEVENT
#define H_PREEVENT 1;

#include "TROOT.h"
// #include "event.hpp"
// #include "TMath.h"

// class preEvent: public event{
//     public:
//         // preEvent();
//         preEvent(Float_t eps_);
//         preEvent(event ev_, Float_t eps_);
//         Float_t Geteps(){return eps;};
//         void   Seteps(Float_t eps_){eps=eps_;};
//         Bool_t DeriveInspect(particle * pT, particle * pP, Bool_t verbose=kFALSE);
//     private:
//         Float_t eps = 1.0 * TMath::DegToRad();
//         Float_t calCor(particle* pT, particle* pP, Bool_t deg=kFALSE);
// };

#endif