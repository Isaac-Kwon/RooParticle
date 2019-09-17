#include "TROOT.h"
// #include "preEvent.hpp"
// #include "path.hpp"
// #include "event.hpp"

// preEvent::preEvent(Float_t eps_):event(){
//     Seteps(eps_);
// }

// preEvent::preEvent(event ev_, Float_t eps_):event(ev_){
//     Seteps(eps_);
// }

// Bool_t preEvent::DeriveInspect(particle * pT, particle * pP,
//                                Bool_t verbose){
//     Int_t i=0;
//     while(kTRUE){
//         if(calCor(pT, pP,kTRUE)<eps){
//             break;
//         }
//         DeriveDT(1);
//         if(verbose && i%1000 == 0){
//             cout<<"i: "<<i<<" | "<< calCor(pT, pP, kFALSE) <<", "<< calCor(pT, pP, kTRUE) <<endl;
//             // pT->GetX().Print();
//             pP->GetX().Print();
//             pP->GetV().Print();

//         }
//         i++;
//     }
//     return kTRUE;  // derivation runs fully
//     // return kFALSE; // derivation doesn't run fully
// }

// Float_t preEvent::calCor(particle* pT, particle* pP, Bool_t deg){
//     TVectorD dX = (pP->GetX())-(pT->GetX());
//     TVectorD dV = (pP->GetV())-(pT->GetV());
//     Double_t cc = (dX*dV)/(TMath::Sqrt(dX.Norm2Sqr() * dV.Norm2Sqr()));
//     if(deg){
//         return TMath::ACos(cc) * TMath::RadToDeg();
//     }
//     return cc;
// }