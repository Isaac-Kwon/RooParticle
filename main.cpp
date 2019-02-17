
#include "TROOT.h"
#include "TVectorT.h"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "eventVolume.hpp"

int main(){

  Int_t i;

  std::cout<<"PARTICLE DEFINIATION"<<endl;

  EMparticle* p1;
  EMparticle* p2;

  Double_t x1_[] = {0., 0., 0.};
  Double_t v1_[] = {0., 0., 0.};

  TVectorD x1 = TVectorD(3, x1_);
  TVectorD v1 = TVectorD(3, v1_);

  Double_t x2_[] = {-10., 0., 0.};
  Double_t v2_[] = {0.04, 0., 0.};

  TVectorD x2 = TVectorD(3, x2_);
  TVectorD v2 = TVectorD(3, v2_);

  p1 = new EMparticle(206,79, x1, v1, true);
  p2 = new EMparticle(1,1, x2, v2);

  // std::cout<<"INITIAL EVENT GENERATION"<<endl;

  event* e = new event();

  // std::cout<<"FORCE DEFINE"<<endl;
  coulombForce* cp = new coulombForce();

  // std::cout<<"EVENT VOLUME GENERATION"<<endl;

  eventVGeneral * vol = new eventVGeneral(40.,40.,40.);

  // std::cout<<"SET PARTICLE"<<endl;

  e->AddParticle(p1);
  e->AddParticle(p2);

  // std::cout<<"SET FORCE"<<endl;

  e->AddForce(cp);

  // std::cout<<"SET VOLUME"<<endl;

  e->SetVolume(vol);

  // std::cout<<"EVENT GENERATION"<<endl;

  TVectorD xx1 = TVectorD(3);
  TVectorD xx2 = TVectorD(3);

  // p1->GetX().Print();
  // p2->GetX().Print();

  for(i=0; i<10; i++){
    std::cout<<"DERIVE NO "<<i<<endl;
    // p1->GetX().Print();
    // p2->GetX().Print();
    xx1 = p1->GetX();
    // xx1.Print();
    xx2 = p2->GetX();
    // xx2.Print();

    std::cout<<"P1 :"<<xx1[0]<<", "<<xx1[1]<<", "<<xx1[2]<<endl;
    std::cout<<"P2 :"<<xx2[0]<<", "<<xx2[1]<<", "<<xx2[2]<<endl;

    e->DeriveDT(1);
  }

  return 0;

}
