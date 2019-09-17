#include "TROOT.h"

// #include "preEvent.hpp"
#include "path.hpp"

#include "event.hpp"
#include "force.hpp"
#include "eventVolume.hpp"
#include "inspector.hpp"

int main(int argc, char** argv){

    //base volume and force
    coulombForce * cp = new coulombForce();
    eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

    Double_t x1_[] = { 0., 0., 0.};
    Double_t x2_[] = { 100., 0., 1000};

    Double_t v1_[] = { 0., 0., 0.};
    Double_t v2_[] = { 0., 0. , -0.05 };

    TVector3 x1 = TVector3(x1_);
    TVector3 x2 = TVector3(x2_);

    TVector3 v1 = TVector3(v1_);
    TVector3 v2 = TVector3(v2_);

    EMparticle * p1 = new EMparticle(196,79, x1, v1, true, false);
    EMparticle * p2 = new EMparticle(4,2, x2, v2, false, true);

    //event template definition
    event * pevent = new event();

    pevent->AddParticle(p1);
    pevent->AddParticle(p2);
    pevent->AddForce(cp);
    pevent->AddInspector(new inspector(p1, p2, "DEG", 0.1));
    pevent->SetVolume(vol);
    pevent->Print(kFALSE, kFALSE, kTRUE);
    pevent->SetInitial();
    pevent->Print(kFALSE, kFALSE, kTRUE);
    pevent->DeriveInspect();
    pevent->Print(kFALSE, kFALSE, kTRUE);


    return 0;
}
