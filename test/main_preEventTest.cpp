#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// #include "preEvent.hpp"
#include "path.hpp"

#include "event.hpp"
#include "force.hpp"
#include "eventVolume.hpp"
#include "inspector.hpp"
#include "TRandom3.h"
class ExperimentT{ //Single Particle Projectile Experiment
public:
  void setTemplate(event* eventT_){eventT=eventT_;}
  virtual void makeEvent(particle* p); // if event already on,
  void offEvent(){eventOn=false;}
  void delEvent(){eventOn=false; delete currentEvent;}
  event* getEvent(){return currentEvent;}
protected:
private:
  event * eventT;
  Bool_t eventOn = false;
  event * currentEvent;
};

void ExperimentT::makeEvent(particle* p){
    if(eventOn){
        this->delEvent();
    }else{
        eventOn = true;
    }
    currentEvent = new event(*eventT);
    currentEvent->AddParticle(p);
}

int main(int argc, char** argv){

    //base volume and force
    coulombForce * cp = new coulombForce();
    eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

    Double_t x1_[] = { 0., 0., 0.};
    TVector3 x1 = TVector3(x1_);
    Double_t v1_[] = { 0., 0., 0.};
    TVector3 v1 = TVector3(v1_);
    
    EMparticle * p01 = new EMparticle(196,79, x1, v1, true, false);
    EMparticle * p02 = new EMparticle(196,79, x1, v1, true, false);

    //event template definition
    event * pevent = new event(); //Normal Event
    pevent->AddParticle(p01);
    pevent->AddForce(cp);
    pevent->SetVolume(vol);

    event * qevent = new event(); //Derive Inspecting Event
    qevent->AddParticle(p02);
    qevent->AddForce(cp);
    qevent->SetVolume(vol);

    ExperimentT * exp1 = new ExperimentT(); //Normal Event (50000 Deriving)
    ExperimentT * exp2 = new ExperimentT(); //Derive Inspecting Experiment
    exp1->setTemplate(pevent);
    exp2->setTemplate(qevent);

    Int_t i;
    
    TRandom3 * r1 = new TRandom3();
    cout<<time(0)<<endl;
    r1->SetSeed(time(0));

    TVector3 x2;
    TVector3 v2;

    TTree * t_ = new TTree("TestTree", "Test for DerivingTime");
    TFile * file_ = new TFile("comparTest.root","RECREATE");
    file_->Add(t_);
    
    Double_t dtime1, dtime2;

    Double_t imp;
    TVector3 x1t, x2t, v1t, v2t;
    t_->Branch("imp", &imp, "imp/D");
    t_->Branch("x1", "TVector3", &x1t);
    t_->Branch("x2", "TVector3", &x2t);
    t_->Branch("v1", "TVector3", &v1t);
    t_->Branch("v2", "TVector3", &v2t);
    t_->Branch("DT1", &dtime1, "DT1/D");
    t_->Branch("DT2", &dtime2, "DT2/D");

    for(i=0; i<5; i++){

        imp = r1->Uniform(0,1000);        
        Double_t x2_[] = { imp, 0., 1000};
        x2 = TVector3(x2_);
        Double_t v2_[] = { 0., 0. , -0.05 };
        v2 = TVector3(v2_);

        EMparticle * p11 = new EMparticle(4,2, x2, v2, false, true);
        EMparticle * p12 = new EMparticle(4,2, x2, v2, false, true);

        cout<<"CR EVENT "<<i<<endl;
        exp1->makeEvent(p11);
        exp2->makeEvent(p12);
        exp1->getEvent()->AddInspector(new inspectorP(exp1->getEvent()->getParticle(0), p11, "DEG", 1.));
        // exp2->getEvent()->AddInspector(new inspector(exp2->getEvent()->getParticle(0), p12, "DEG", 0.1));

        // exp2->getEvent()->Print();
        exp1->getEvent()->DeriveInspect(1);
        exp2->getEvent()->DeriveDTN(50000);

        x1t = exp1->getEvent()->getParticle(1)->GetX();
        x2t = exp2->getEvent()->getParticle(1)->GetX();
        v1t = exp1->getEvent()->getParticle(1)->GetV();
        v2t = exp2->getEvent()->getParticle(1)->GetV();
        dtime1 = exp1->getEvent()->getParticle(1)->GetTime();
        dtime2 = exp2->getEvent()->getParticle(1)->GetTime();
        
        t_->Fill();
        t_->AutoSave();

        exp1->delEvent();
        exp2->delEvent();

    }


    file_->Write();
    file_->Close();
    


    return 0;
}
