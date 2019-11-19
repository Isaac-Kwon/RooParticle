#include "experiment.hpp"
#include "particle.hpp"
#include "event.hpp"

#include "TROOT.h"

#include "cstdarg"

void experiment::MakeEvent(Int_t ninparticle, ...){
    eventOn = kTRUE;
    currentEvent = new event(*eventT);
    Int_t i;
    particle* p;
    va_list argList;
    va_start(argList, ninparticle);

    for(i=0; i<ninparticle; i++){
        p = va_arg(argList, particle* );
        currentEvent->AddParticle(p);
    }
    va_end(argList);

    return;
}

void experimentP::MakeEvent(particle* p1, ...){
    eventOn = kTRUE;
    currentEvent = new event(*eventT);

    particle* p=p1;
    va_list argList;
    currentEvent->AddParticle(p1);
    va_start(argList, p1);

    Int_t i;
    for(i=0; i<(ninparticle-1); i++){
        p = va_arg(argList, particle* );
        currentEvent->AddParticle(p);
    }
    va_end(argList);

    return;
}
