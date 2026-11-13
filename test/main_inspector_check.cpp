////////////////////////////////////////////
//
// RooParticle Test Application
// Inspector check (inspectorP and inspectorE)
//
// Every method is checked against numbers calculated by hand:
// the Evaluate() value, and the direction (strict > or <) of Inspect().
// Also checked : angles of 1e-9 rad and (pi - 1e-9) rad (where acos fails),
// zero dX / dV, string <-> enum construction and unknown method names.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "eventVolume.hpp"
#include "inspector.hpp"

#include "iostream"
#include <cstdio>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

static Int_t gNCheck = 0;
static Int_t gNFail  = 0;

void Check(const char * what, Bool_t ok){
  gNCheck++;
  if(!ok){
    gNFail++;
    printf("  FAIL  %s\n", what);
  }
}

// |got - want| <= atol + rtol |want|
void CheckNear(const char * what, Double_t got, Double_t want, Double_t rtol, Double_t atol = 0.){
  gNCheck++;
  if(!(TMath::Abs(got-want) <= atol + rtol*TMath::Abs(want))){
    gNFail++;
    printf("  FAIL  %s : got %.17g, want %.17g\n", what, got, want);
  }
}

// A force with fixed numbers : acceleration 5/m (3,4,0 direction) and potential -2 per pair.
class fixedForce: public force{
public:
  fixedForce():force(){}
  TVector3 Evaluate(particle* p1_, particle* p2_){return TVector3(3.,4.,0.);}
  Double_t Potential(particle* p1_, particle* p2_){return -2.;}
};

// =============================== inspectorP ===============================

// p1 - p2 = dX , v1 - v2 = dV
void SetPair(particle & p1, particle & p2, TVector3 dX, TVector3 dV){
  p1.applyMechanic(dX, dV);
  p2.applyMechanic(TVector3(0.,0.,0.), TVector3(0.,0.,0.));
}

// Evaluate() and the direction of Inspect() around the value :
// greater=kTRUE  : satisfied when Evaluate() > val
// greater=kFALSE : satisfied when Evaluate() < val
void CheckP(particle & p1, particle & p2, const char * name, Double_t expect,
            Bool_t greater, Double_t margin, Double_t rtol = 1e-12){
  char buf[200];
  inspectorP * ins;

  ins = new inspectorP(&p1, &p2, name, expect);
  snprintf(buf, 200, "P %s Evaluate", name);
  CheckNear(buf, ins->Evaluate(), expect, rtol, 1e-300);
  // equal to the evaluated value : strict comparison, never satisfied
  ins->SetValue(ins->Evaluate());
  snprintf(buf, 200, "P %s Inspect at the value (strict)", name);
  Check(buf, !ins->Inspect());
  delete ins;

  ins = new inspectorP(&p1, &p2, name, expect - margin);
  snprintf(buf, 200, "P %s Inspect, value below", name);
  Check(buf, ins->Inspect() == greater);
  delete ins;

  ins = new inspectorP(&p1, &p2, name, expect + margin);
  snprintf(buf, 200, "P %s Inspect, value above", name);
  Check(buf, ins->Inspect() == !greater);
  delete ins;
}

void CheckPairMethods(){
  printf("== inspectorP : angle, cosine ==\n");
  particle p1(1., TVector3(), TVector3()), p2(1., TVector3(), TVector3());

  // dX = (1,0,0) ; dV, expected angle (rad)
  struct {TVector3 dV; Double_t rad; const char * what;} cases[] = {
    {TVector3(1.,0.,0.),   0.,           "0 degree"},
    {TVector3(0.,2.,0.),   TMath::Pi()/2., "90 degree"},
    {TVector3(-3.,0.,0.),  TMath::Pi(),   "180 degree"},
    {TVector3(1.,1.,0.),   TMath::Pi()/4., "45 degree"},
    {TVector3(-1.,1.,0.),  3.*TMath::Pi()/4., "135 degree"},
    {TVector3(1.,-1.,0.),  TMath::Pi()/4., "45 degree (other side)"},
    {TVector3(0.,0.,5.),   TMath::Pi()/2., "90 degree (z)"},
  };
  for(auto & c : cases){
    printf("  dX=(1,0,0) dV=(%g,%g,%g) : %s\n", c.dV.X(), c.dV.Y(), c.dV.Z(), c.what);
    SetPair(p1, p2, TVector3(1.,0.,0.), c.dV);
    char buf[100];
    inspectorP rad(&p1, &p2, "RAD", 0.), deg(&p1, &p2, "DEG", 0.), cs(&p1, &p2, "COS", 0.);
    snprintf(buf, 100, "RAD %s", c.what);
    CheckNear(buf, rad.Evaluate(), c.rad, 1e-14, 1e-15);
    snprintf(buf, 100, "DEG %s", c.what);
    CheckNear(buf, deg.Evaluate(), c.rad*TMath::RadToDeg(), 1e-14, 1e-13);
    snprintf(buf, 100, "COS %s", c.what);
    CheckNear(buf, cs.Evaluate(), TMath::Cos(c.rad), 1e-14, 1e-15);
  }

  // directions : 45 degree
  SetPair(p1, p2, TVector3(1.,0.,0.), TVector3(1.,1.,0.));
  CheckP(p1, p2, "DEG", 45., kFALSE, 0.5);
  CheckP(p1, p2, "RAD", TMath::Pi()/4., kFALSE, 0.01);
  CheckP(p1, p2, "COS", TMath::Sqrt(0.5), kFALSE, 0.01, 1e-14);

  // the precision cases
  printf("== inspectorP : 1e-9 rad and (pi - 1e-9) rad ==\n");
  const Double_t eps = 1e-9;
  SetPair(p1, p2, TVector3(1.,0.,0.), TVector3(1.,eps,0.));
  {
    inspectorP rad(&p1, &p2, "RAD", 0.), deg(&p1, &p2, "DEG", 0.), cs(&p1, &p2, "COS", 0.);
    printf("  small : RAD = %.12e (acos would give %.3e)\n", rad.Evaluate(), TMath::ACos(cs.Evaluate()));
    CheckNear("RAD 1e-9", rad.Evaluate(), eps, 1e-12);
    CheckNear("DEG 1e-9 rad", deg.Evaluate(), eps*TMath::RadToDeg(), 1e-12);
    Check("COS 1e-9 rad is 1", cs.Evaluate() <= 1. && cs.Evaluate() >= 1.-1e-15);
    Check("acos(cos) lost the angle (why atan2)", TMath::ACos(cs.Evaluate()) == 0.);
    inspectorP a(&p1, &p2, "RAD", 2e-9), b(&p1, &p2, "RAD", 5e-10), c(&p1, &p2, "DEG", 1.);
    Check("RAD < 2e-9 for 1e-9", a.Inspect());
    Check("RAD < 5e-10 not for 1e-9", !b.Inspect());
    Check("DEG < 1 for 1e-9 rad", c.Inspect());
  }
  SetPair(p1, p2, TVector3(1.,0.,0.), TVector3(-1.,eps,0.));
  {
    inspectorP rad(&p1, &p2, "RAD", 0.), deg(&p1, &p2, "DEG", 0.), cs(&p1, &p2, "COS", 0.);
    Double_t r = rad.Evaluate();
    printf("  large : pi - RAD = %.12e (acos would give pi - %.3e)\n", TMath::Pi()-r, TMath::Pi()-TMath::ACos(cs.Evaluate()));
    CheckNear("RAD pi-1e-9", r, TMath::Pi()-eps, 0., 1e-15);
    CheckNear("pi - RAD", TMath::Pi()-r, eps, 1e-6);
    CheckNear("DEG pi-1e-9", deg.Evaluate(), (TMath::Pi()-eps)*TMath::RadToDeg(), 1e-14);
    Check("COS pi-1e-9 rad is -1", cs.Evaluate() >= -1. && cs.Evaluate() <= -1.+1e-15);
    Check("acos(cos) lost the angle (why atan2)", TMath::ACos(cs.Evaluate()) == TMath::Pi());
    inspectorP a(&p1, &p2, "DEG", 179.99999999), b(&p1, &p2, "DEG", 179.);
    // 179.99999999 degree = pi - 1.7e-10 rad is above pi - 1e-9 : satisfied
    Check("DEG < 179.99999999 for pi-1e-9", a.Inspect());
    Check("DEG < 179 not for pi-1e-9", !b.Inspect());
  }

  // cosine never leaves [-1,1], also with rounding (parallel and anti-parallel vectors)
  printf("== inspectorP : rounding of the cosine ==\n");
  TRandom3 rnd(12345);
  Int_t nbad = 0, nrange = 0;
  for(Int_t i=0; i<2000; i++){
    TVector3 dX(rnd.Uniform(-10,10), rnd.Uniform(-10,10), rnd.Uniform(-10,10));
    Double_t k = rnd.Uniform(0.1,10.) * (i%2 ? 1. : -1.);
    SetPair(p1, p2, dX, dX*k);
    inspectorP cs(&p1, &p2, "COS", 0.), rad(&p1, &p2, "RAD", 0.);
    Double_t c = cs.Evaluate(), r = rad.Evaluate();
    if(!(c>=-1. && c<=1.)) nrange++;
    Double_t expect = (k>0) ? 0. : TMath::Pi();
    // atan2 of the (rounded) cross product : 1e-15 level
    if(!(TMath::Abs(r-expect) < 1e-14)) nbad++;
  }
  Check("cosine of parallel / anti-parallel vectors stays in [-1,1]", nrange==0);
  Check("angle of parallel / anti-parallel vectors is 0 / pi (< 1e-14)", nbad==0);

  printf("== inspectorP : speed, distance, radial velocity, count ==\n");
  // dX = (3,4,0) : |dX| = 5 ; dV = (1,2,0) : |dV| = sqrt5, dX.dV = 11, VR = 2.2
  SetPair(p1, p2, TVector3(3.,4.,0.), TVector3(1.,2.,0.));
  CheckP(p1, p2, "DSTG", 5., kTRUE,  0.5);
  CheckP(p1, p2, "DSTL", 5., kFALSE, 0.5);
  CheckP(p1, p2, "SPDG", TMath::Sqrt(5.), kTRUE,  0.5);
  CheckP(p1, p2, "SPDL", TMath::Sqrt(5.), kFALSE, 0.5);
  CheckP(p1, p2, "VRG",  2.2, kTRUE,  0.5);
  CheckP(p1, p2, "VRL",  2.2, kFALSE, 0.5);
  {
    // radial velocity, sign : approaching is negative
    inspectorP g(&p1, &p2, "VRG", 0.), l(&p1, &p2, "VRL", 0.);
    Check("receding : VRG>0", g.Inspect());
    Check("receding : VRL<0 false", !l.Inspect());
    SetPair(p1, p2, TVector3(3.,4.,0.), TVector3(-1.,-2.,0.));
    CheckNear("VR approaching", g.Evaluate(), -2.2, 1e-14);
    Check("approaching : VRG>0 false", !g.Inspect());
    Check("approaching : VRL<0", l.Inspect());
  }
  // the particle order : dX, dV are p1 - p2
  SetPair(p1, p2, TVector3(0.,3.,4.), TVector3(0.,0.,0.));
  {
    inspectorP a(&p1, &p2, "DSTL", 6.), b(&p2, &p1, "DSTL", 6.);
    CheckNear("DST from swapped particles", b.Evaluate(), 5., 1e-14);
    Check("DSTL<6 both orders", a.Inspect() && b.Inspect());
  }
  // CNT : the smaller step count of the two (here 3 and 5)
  {
    for(Int_t i=0; i<5; i++) p1.finishStep(1.);
    for(Int_t i=0; i<3; i++) p2.finishStep(1.);
    CheckP(p1, p2, "CNT", 3., kTRUE, 1.);
  }

  printf("== inspectorP : zero dX or dV (never satisfied) ==\n");
  particle q1(1., TVector3(), TVector3()), q2(1., TVector3(), TVector3());
  const char * angular[] = {"COS", "DEG", "RAD"};
  // dV = 0 : angle undefined, even with a criterion that any angle would pass
  SetPair(q1, q2, TVector3(1.,2.,3.), TVector3(0.,0.,0.));
  for(const char * n : angular){
    inspectorP ins(&q1, &q2, n, 1e9);
    char buf[100];
    snprintf(buf, 100, "dV=0 : %s Evaluate is NaN", n);
    Check(buf, TMath::IsNaN(ins.Evaluate()));
    snprintf(buf, 100, "dV=0 : %s Inspect is false", n);
    Check(buf, !ins.Inspect());
  }
  // dX = 0
  SetPair(q1, q2, TVector3(0.,0.,0.), TVector3(1.,2.,3.));
  for(const char * n : angular){
    inspectorP ins(&q1, &q2, n, 1e9);
    char buf[100];
    snprintf(buf, 100, "dX=0 : %s Evaluate is NaN", n);
    Check(buf, TMath::IsNaN(ins.Evaluate()));
    snprintf(buf, 100, "dX=0 : %s Inspect is false", n);
    Check(buf, !ins.Inspect());
  }
  {
    inspectorP g(&q1, &q2, "VRG", -1e9), l(&q1, &q2, "VRL", 1e9);
    Check("dX=0 : VRG Evaluate is NaN", TMath::IsNaN(g.Evaluate()));
    Check("dX=0 : VRG, VRL Inspect are false", !g.Inspect() && !l.Inspect());
    // methods that do not need a direction still work
    inspectorP d(&q1, &q2, "DSTL", 1.), s(&q1, &q2, "SPDG", 1.);
    Check("dX=0 : DSTL<1 is satisfied (distance 0)", d.Inspect());
    Check("dX=0 : SPDG>1 is satisfied", s.Inspect());
  }
  // both zero
  SetPair(q1, q2, TVector3(0.,0.,0.), TVector3(0.,0.,0.));
  {
    inspectorP deg(&q1, &q2, "DEG", 1.);
    Check("dX=dV=0 : DEG Inspect is false (atan2(0,0) must not mean 0 degree)", !deg.Inspect());
  }
}

// =============================== inspectorE ===============================

// Value, G side and L side of an inspectorE method, around the expected value.
void CheckE(event & ev, const char * name, Double_t expect, Bool_t greater, Double_t margin){
  char buf[200];
  inspectorE * ins;

  ins = new inspectorE(&ev, name, expect);
  snprintf(buf, 200, "E %s Evaluate", name);
  CheckNear(buf, ins->Evaluate(), expect, 1e-12, 1e-300);
  snprintf(buf, 200, "E %s Inspect at the value (strict)", name);
  Check(buf, !ins->Inspect());
  delete ins;

  ins = new inspectorE(&ev, name, expect - margin);
  snprintf(buf, 200, "E %s Inspect, value below", name);
  Check(buf, ins->Inspect() == greater);
  delete ins;

  ins = new inspectorE(&ev, name, expect + margin);
  snprintf(buf, 200, "E %s Inspect, value above", name);
  Check(buf, ins->Inspect() == !greater);
  delete ins;
}

void CheckEventMethods(){
  printf("== inspectorE ==\n");
  // m1=2 at origin with v=(3,0,0), m2=4 at (3,4,0) at rest ; forces and potential fixed (fixedForce)
  //   KE = 0.5*2*9 = 9      PE = -2 (pair, each particle sees -2, half-counted)
  //   |PE/KE| = 2/9         E = 7
  //   min distance 5        a1 = 5/2 = 2.5 , a2 = 5/4 = 1.25
  //   steps (p1,p2) = (5,2) time (p1,p2) = (5,2)
  fixedForce ff;
  event ev;
  ev.AddForce(&ff);
  particle * p1 = ev.AddParticle(new particle(2., TVector3(0.,0.,0.), TVector3(3.,0.,0.), false, false));
  particle * p2 = ev.AddParticle(new particle(4., TVector3(3.,4.,0.), TVector3(0.,0.,0.), false, false));
  for(Int_t i=0; i<5; i++) p1->finishStep(1.);
  for(Int_t i=0; i<2; i++) p2->finishStep(1.);
  eventVGeneral box(-10.,10.,-10.,10.,-10.,10.);
  ev.SetVolume(&box);

  CheckNear("event net KE (setup)", ev.GetNetKE(), 9., 1e-14);
  CheckNear("event net PE (setup)", ev.GetNetPE(), -2., 1e-14);

  CheckE(ev, "KEG", 9., kTRUE,  1.);
  CheckE(ev, "KEL", 9., kFALSE, 1.);
  CheckE(ev, "PEG", -2., kTRUE,  1.);
  CheckE(ev, "PEL", -2., kFALSE, 1.);
  CheckE(ev, "PKG", 2./9., kTRUE,  0.1);
  CheckE(ev, "PKL", 2./9., kFALSE, 0.1);
  CheckE(ev, "TMG", 5., kTRUE,  1.);
  CheckE(ev, "TML", 5., kFALSE, 1.);
  CheckE(ev, "DSG", 5., kTRUE,  1.);
  CheckE(ev, "DSL", 5., kFALSE, 1.);
  CheckE(ev, "ACG", 2.5, kTRUE,  0.5);
  CheckE(ev, "ACL", 2.5, kFALSE, 0.5);

  // step count : Evaluate DX* = largest (5), DN* = smallest (2).
  // Inspect (existing rule) : DXG, DNG compare the largest ; DXL, DNL the smallest.
  {
    const char * names[] = {"DXG", "DXL", "DNG", "DNL"};
    Double_t evalue[] = {5., 5., 2., 2.};
    for(Int_t i=0; i<4; i++){
      char buf[100];
      inspectorE ins(&ev, names[i], 0.);
      snprintf(buf, 100, "E %s Evaluate", names[i]);
      CheckNear(buf, ins.Evaluate(), evalue[i], 0.);
    }
    inspectorE a(&ev,"DXG",4.), b(&ev,"DXG",5.), c(&ev,"DNG",4.), d(&ev,"DNG",5.);
    Check("DXG < largest(5): val 4 satisfied, val 5 not", a.Inspect() && !b.Inspect());
    Check("DNG < largest(5): val 4 satisfied, val 5 not", c.Inspect() && !d.Inspect());
    inspectorE e(&ev,"DXL",3.), f(&ev,"DXL",2.), g(&ev,"DNL",3.), h(&ev,"DNL",2.);
    Check("DXL smallest(2) < val: val 3 satisfied, val 2 not", e.Inspect() && !f.Inspect());
    Check("DNL smallest(2) < val: val 3 satisfied, val 2 not", g.Inspect() && !h.Inspect());
  }

  // energy drift : E0 = 7 is taken at the first Evaluate()
  {
    inspectorE g(&ev, "EDG", 0.5), l(&ev, "EDL", 0.5);
    CheckNear("EDG at the first inspection is 0", g.Evaluate(), 0., 0., 1e-300);
    Check("EDG>0.5 not yet / EDL<0.5 yes", !g.Inspect() && l.Inspect());
    p1->applyV(TVector3(4.,0.,0.));      // KE = 16 , E = 14 : drift |14-7|/7 = 1
    CheckNear("EDG after the energy doubled", g.Evaluate(), 1., 1e-14);
    Check("EDG>0.5 yes / EDL<0.5 not", g.Inspect() && !l.Inspect());
    inspectorE g2(&ev, "EDG", 1.);
    CheckNear("fresh EDG takes the reference now (E0=14)", g2.Evaluate(), 0., 0., 1e-300);
    g.ResetReference();
    CheckNear("ResetReference", g.Evaluate(), 0., 0., 1e-300);
    p1->applyV(TVector3(3.,0.,0.));
  }

  // PK with zero kinetic energy : the largest number, "G" satisfied
  {
    particle * q1 = ev.getParticle(0);
    TVector3 v1 = q1->GetV();
    q1->applyV(TVector3());
    inspectorE g(&ev, "PKG", 1e300), l(&ev, "PKL", 1e300);
    Check("PK with KE=0 : PKG satisfied, PKL not", g.Inspect() && !l.Inspect());
    q1->applyV(v1);
  }

  // volume
  printf("== inspectorE : volume ==\n");
  {
    inspectorE vo(&ev, "VOUT"), va(&ev, "VOUTA");
    Check("both inside : VOUT false", !vo.Inspect() && vo.Evaluate()==0.);
    Check("both inside : VOUTA false", !va.Inspect() && va.Evaluate()==0.);
    p2->applyX(TVector3(50.,0.,0.));
    Check("one out : VOUT true", vo.Inspect() && vo.Evaluate()==1.);
    Check("one out : VOUTA false", !va.Inspect() && va.Evaluate()==0.);
    p1->applyX(TVector3(0.,50.,0.));
    Check("all out : VOUT true", vo.Inspect());
    Check("all out : VOUTA true", va.Inspect() && va.Evaluate()==1.);
  }
  // invincible particles are ignored by the volume
  {
    event evb;
    particle * m = evb.AddParticle(new particle(1., TVector3(0.,0.,0.), TVector3(), false, false));
    evb.AddParticle(new particle(1., TVector3(1000.,0.,0.), TVector3(), true, false));
    eventVGeneral box2(-10.,10.,-10.,10.,-10.,10.);
    evb.SetVolume(&box2);
    inspectorE vo(&evb, "VOUT"), va(&evb, "VOUTA");
    Check("invincible outside is ignored : VOUT false, VOUTA false", !vo.Inspect() && !va.Inspect());
    m->applyX(TVector3(100.,0.,0.));
    Check("moving particle out : VOUT true, VOUTA true", vo.Inspect() && va.Inspect());
  }
}

// ==================== string / enum construction, invalid names ====================

Bool_t SameValue(Double_t a, Double_t b){
  return (TMath::IsNaN(a) && TMath::IsNaN(b)) || a==b;
}

void CheckConstruction(){
  printf("== string and enum construction ==\n");
  particle p1(2., TVector3(3.,4.,1.), TVector3(1.,2.,0.5)), p2(1., TVector3(0.,0.,0.), TVector3(-1.,0.,0.));
  for(Int_t i=0; i<4; i++) p1.finishStep(1.);

  const char * pn[] = {"COS","DEG","RAD","CNT","SPDG","SPDL","DSTG","DSTL","VRG","VRL"};
  for(Int_t i=0; i<10; i++){
    inspector::methodID id = (inspector::methodID)(inspector::kCOS + i);
    inspectorP s(&p1, &p2, pn[i], 1.5), e(&p1, &p2, id, 1.5);
    char buf[100];
    snprintf(buf, 100, "P %s : enum == string (value, result, name, validity)", pn[i]);
    Check(buf, SameValue(s.Evaluate(), e.Evaluate()) && s.Inspect()==e.Inspect()
            && s.GetMethod()==e.GetMethod() && e.GetMethod()==pn[i]
            && s.GetMethodID()==id && e.GetMethodID()==id
            && s.IsValid() && e.IsValid());
    Check("P is inspectorP only", s.IsInspectorP() && !s.IsInspectorE());
    // SetMethod changes both
    inspectorP t(&p1, &p2, "DEG", 1.5);
    t.SetMethod(id);
    Check("P SetMethod(enum) follows", t.GetMethod()==pn[i] && SameValue(t.Evaluate(), s.Evaluate()));
    t.SetMethod("DEG");
    t.SetMethod(TString(pn[i]));
    Check("P SetMethod(TString) follows", t.GetMethodID()==id && SameValue(t.Evaluate(), s.Evaluate()));
  }

  event ev;
  fixedForce ff;
  ev.AddForce(&ff);
  ev.AddParticle(new particle(2., TVector3(0.,0.,0.), TVector3(3.,0.,0.)));
  ev.AddParticle(new particle(4., TVector3(3.,4.,0.), TVector3(0.,0.,0.)));
  eventVGeneral box(-10.,10.,-10.,10.,-10.,10.);
  ev.SetVolume(&box);
  const char * en[] = {"KEG","KEL","PEG","PEL","PKG","PKL","EDG","EDL","TMG","TML","DXG","DXL","DNG","DNL",
                       "DSG","DSL","ACG","ACL","VOUT","VOUTA"};
  for(Int_t i=0; i<20; i++){
    inspector::methodID id = (inspector::methodID)(inspector::kKEG + i);
    inspectorE s(&ev, en[i], 3.), e(&ev, id, 3.);
    char buf[100];
    snprintf(buf, 100, "E %s : enum == string (value, result, name, validity)", en[i]);
    Check(buf, SameValue(s.Evaluate(), e.Evaluate()) && s.Inspect()==e.Inspect()
            && s.GetMethod()==e.GetMethod() && e.GetMethod()==en[i]
            && s.GetMethodID()==id && e.GetMethodID()==id
            && s.IsValid() && e.IsValid());
    Check("E is inspectorE only", s.IsInspectorE() && !s.IsInspectorP());
  }
  Check("enum value is kept", inspector::MethodName(inspector::kVOUTA)=="VOUTA"
        && inspector::ParseMethod("VOUTA")==inspector::kVOUTA
        && inspector::MethodName(inspector::kInvalid)=="");

  printf("== unknown names ==\n");
  const char * bad[] = {"DGE", "", "deg", "DEG ", "KEG", "VOUTX", "COS2"};   // "KEG" is not a pair method
  for(const char * b : bad){
    inspectorP ip(&p1, &p2, b, 1.5);
    char buf[100];
    snprintf(buf, 100, "inspectorP \"%s\" invalid", b);
    Check(buf, !ip.IsValid() && !ip.Inspect() && ip.GetMethod()==b);
  }
  const char * badE[] = {"DGE", "", "keg", "DEG", "COS", "VOUTAA"};   // "DEG" is not an event method
  for(const char * b : badE){
    inspectorE ie(&ev, b, 1.5);
    char buf[100];
    snprintf(buf, 100, "inspectorE \"%s\" invalid", b);
    Check(buf, !ie.IsValid() && !ie.Inspect() && ie.GetMethod()==b);
  }
  {
    inspectorP ip(&p1, &p2, inspector::kKEG, 1.5);     // enum of the other kind
    inspectorE ie(&ev, inspector::kDEG, 1.5);
    inspectorP in(&p1, &p2, inspector::kInvalid, 1.5);
    Check("enum of the other kind / kInvalid is invalid", !ip.IsValid() && !ie.IsValid() && !in.IsValid()
          && !ip.Inspect() && !ie.Inspect() && !in.Inspect());
  }
  {
    // a valid inspector made invalid and valid again by SetMethod
    inspectorP ip(&p1, &p2, "DSTG", 1.);
    Check("valid", ip.IsValid());
    ip.SetMethod("FOO");
    Check("renamed to FOO : invalid, name kept", !ip.IsValid() && ip.GetMethod()=="FOO" && !ip.Inspect());
    ip.SetMethod("DSTL");
    Check("renamed to DSTL : valid again", ip.IsValid() && ip.GetMethod()=="DSTL");
  }

  printf("== common interface ==\n");
  {
    inspectorP ip(&p1, &p2, "DEG", 1.);
    Check("criterion", ip.GetCriterion()==1.);
    ip.SetValue(2.5);
    Check("SetValue", ip.GetCriterion()==2.5);
    Check("abort is off by default", !ip.IsAbort());
    ip.SetAbort();
    Check("SetAbort()", ip.IsAbort());
    ip.SetAbort(kFALSE);
    Check("SetAbort(kFALSE)", !ip.IsAbort());
    inspector * base = &ip;
    Check("virtual IsInspectorP/E", base->IsInspectorP() && !base->IsInspectorE());
    // as the event copies it : through the base class
    inspectorP copy(&p1, &p2, base->GetMethod(), base->GetCriterion());
    Check("copy through GetMethod()/GetCriterion()", copy.GetMethodID()==inspector::kDEG && copy.GetCriterion()==2.5);
  }
}

int main(){
  CheckPairMethods();
  CheckEventMethods();
  CheckConstruction();

  printf("\n%d checks, %d failed\n", gNCheck, gNFail);
  if(gNFail){
    printf("INSPECTOR CHECK FAILED\n");
    return 1;
  }
  printf("ALL CHECKS PASSED\n");
  return 0;
}
