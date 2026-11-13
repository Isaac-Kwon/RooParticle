////////////////////////////////////////////
//
// RooParticle Test Application
// Force and particle check (no argument)
//
// 1. Every force class : force = -grad(potential) by central finite
//    differences at several separations, action = -reaction, sign.
// 2. Softening : 0 is bit-identical to no softening (single evaluation and
//    a whole trajectory), > 0 stays finite at r -> 0.
// 3. Record period : number of path points and recorded values.
// 4. Particle copy : step count, charge, record period.
// 5. Small helpers : GetGamma, GetMomentum, SetK / SetG.
//
// Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "RPConfig.hpp"

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"

#include "iostream"
#include <cstdio>
#include <cstring>
#include <vector>

#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"

static Int_t gFail = 0;

// One row of the table : value compared with a tolerance (value <= tol passes).
void Row(const TString name, Double_t value, Double_t tol){
  Bool_t ok = TMath::Finite(value) && value<=tol;
  if(!ok) gFail++;
  printf("  %-58s %12.3e %10.1e  %s\n", name.Data(), value, tol, ok ? "PASS" : "FAIL");
}

// One row for a yes / no check.
void RowB(const TString name, Bool_t ok, const TString detail = ""){
  if(!ok) gFail++;
  printf("  %-58s %23s  %s\n", name.Data(), detail.Data(), ok ? "PASS" : "FAIL");
}

void Header(const char* title){
  printf("\n%s\n", title);
  printf("  %-58s %12s %10s  %s\n", "check", "value", "tol", "result");
}

Bool_t SameBits(Double_t a, Double_t b){return memcmp(&a, &b, sizeof(Double_t))==0;}
Bool_t SameBits(TVector3 a, TVector3 b){return SameBits(a.X(),b.X()) && SameBits(a.Y(),b.Y()) && SameBits(a.Z(),b.Z());}

struct forceCase{
  TString  name;
  force*   f;
  Bool_t   attractive;    // expected sign for the two test particles (like charges)
  std::vector<Double_t> r; // separations [fm]
};

// Unit direction that is not on an axis, so that every component is tested.
const TVector3 gDir = TVector3(0.6, -0.48, 0.64).Unit();
const TVector3 gX2  = TVector3(1.3, -0.7, 2.1);

// Largest relative error between Evaluate() and -grad Potential() by central
// differences on the position of p1, for p1 at distance r from p2.
Double_t GradientError(force* f, particle* p1, particle* p2, Double_t r){
  TVector3 x1 = gX2 + r*gDir;
  p1->applyX(x1);
  TVector3 F = f->Evaluate(p1, p2);
  Double_t h = 1e-4 * r;
  TVector3 Fn;
  Int_t i;
  for(i=0; i<3; i++){
    TVector3 e; e(i) = h;
    p1->applyX(x1 + e);
    Double_t up = f->Potential(p1, p2);
    p1->applyX(x1 - e);
    Double_t um = f->Potential(p1, p2);
    Fn(i) = -(up - um)/(2*h);
  }
  p1->applyX(x1);
  return (F - Fn).Mag() / F.Mag();
}

void CheckForces(){
  coulombForce    coul;
  coulombForce    coulS;  coulS.SetSoftening(2.);
  universalForce  grav;
  universalForce  gravS;  gravS.SetSoftening(2.);
  yukawaForce     yuk;
  woodsSaxonForce ws;

  std::vector<forceCase> cases;
  Double_t rC[]  = {0.5, 1., 3., 10., 100., 1000.};
  Double_t rY[]  = {0.3, 1., 1.4138, 3., 6., 12.};
  Double_t rW[]  = {1., 5., 8., 8.886, 10., 15.};
  forceCase c;
  c.name="coulomb";             c.f=&coul;  c.attractive=kFALSE; c.r.assign(rC, rC+6); cases.push_back(c);
  c.name="coulomb eps=2";       c.f=&coulS; c.attractive=kFALSE; c.r.assign(rC, rC+6); cases.push_back(c);
  c.name="gravity";             c.f=&grav;  c.attractive=kTRUE;  c.r.assign(rC, rC+6); cases.push_back(c);
  c.name="gravity eps=2";       c.f=&gravS; c.attractive=kTRUE;  c.r.assign(rC, rC+6); cases.push_back(c);
  c.name="yukawa";              c.f=&yuk;   c.attractive=kTRUE;  c.r.assign(rY, rY+6); cases.push_back(c);
  c.name="woods-saxon";         c.f=&ws;    c.attractive=kTRUE;  c.r.assign(rW, rW+6); cases.push_back(c);

  // alpha and Au, like charges
  EMparticle pa(4, 2, TVector3(), TVector3(), false, false);
  EMparticle pb(197, 79, gX2, TVector3(), false, false);

  Header("1. force = -grad(potential), relative error (central differences, h = 1e-4 r)");
  size_t k, j;
  for(k=0; k<cases.size(); k++){
    for(j=0; j<cases[k].r.size(); j++){
      Double_t r = cases[k].r[j];
      Row(TString::Format("%-14s r = %8.4f fm", cases[k].name.Data(), r),
          GradientError(cases[k].f, &pa, &pb, r), 1e-6);
    }
  }

  Header("2. action = -reaction : |F(a on b) + F(b on a)| / |F|, |U(a,b) - U(b,a)| / |U|");
  for(k=0; k<cases.size(); k++){
    Double_t worstF = 0., worstU = 0.;
    for(j=0; j<cases[k].r.size(); j++){
      pa.applyX(gX2 + cases[k].r[j]*gDir);
      TVector3 Fab = cases[k].f->Evaluate(&pb, &pa);
      TVector3 Fba = cases[k].f->Evaluate(&pa, &pb);
      Double_t Uab = cases[k].f->Potential(&pb, &pa);
      Double_t Uba = cases[k].f->Potential(&pa, &pb);
      worstF = TMath::Max(worstF, (Fab + Fba).Mag() / Fab.Mag());
      worstU = TMath::Max(worstU, TMath::Abs(Uab - Uba) / TMath::Abs(Uab));
    }
    Row(TString::Format("%-14s force", cases[k].name.Data()), worstF, 1e-14);
    Row(TString::Format("%-14s potential", cases[k].name.Data()), worstU, 1e-14);
  }

  Header("3. sign of the force (radial component F.(x1-x2)) and of the potential");
  for(k=0; k<cases.size(); k++){
    Bool_t ok = kTRUE;
    for(j=0; j<cases[k].r.size(); j++){
      pa.applyX(gX2 + cases[k].r[j]*gDir);
      Double_t fr = cases[k].f->Evaluate(&pa, &pb).Dot(pa.GetX() - pb.GetX());
      Double_t u  = cases[k].f->Potential(&pa, &pb);
      if(cases[k].attractive) ok = ok && fr<0 && u<0;
      else                    ok = ok && fr>0 && u>0;
    }
    RowB(TString::Format("%-14s %s at every r", cases[k].name.Data(),
         cases[k].attractive ? "attracts" : "repels (like charges)"), ok);
  }
  { // opposite charges attract
    EMparticle pe(1, -1, gX2 + 10.*gDir, TVector3(), false, false);
    Double_t fr = coul.Evaluate(&pe, &pb).Dot(pe.GetX() - pb.GetX());
    RowB("coulomb        opposite charges attract", fr<0 && coul.Potential(&pe, &pb)<0);
  }
  { // a neutral particle feels no Coulomb force (charge read from particle*)
    particle pn(1, gX2 + 10.*gDir);
    RowB("coulomb        neutral particle : no force", coul.Evaluate(&pn, &pb).Mag()==0.);
  }
  { // negative Yukawa coupling is repulsive
    yukawaForce yr(-0.1);
    pa.applyX(gX2 + 2.*gDir);
    RowB("yukawa g < 0   repels", yr.Evaluate(&pa, &pb).Dot(pa.GetX() - pb.GetX())>0);
  }
  { // woods-saxon : finite everywhere, half depth at r = R
    pa.applyX(gX2 + ws.GetRadius()*gDir);
    Row("woods-saxon    |U(R) + V0/2| / V0", TMath::Abs(ws.Potential(&pa, &pb) + 0.5*ws.GetDepth())/ws.GetDepth(), 1e-15);
    pa.applyX(gX2);
    RowB("woods-saxon    finite at r = 0",
         TMath::Finite(ws.Potential(&pa, &pb)) && ws.Evaluate(&pa, &pb).Mag()==0.);
    pa.applyX(gX2 + 1e4*gDir);
    RowB("woods-saxon    finite at r = 1e4 fm (no overflow)",
         TMath::Finite(ws.Potential(&pa, &pb)) && TMath::Finite(ws.Evaluate(&pa, &pb).Mag()));
  }

  Header("4. parameters");
  RowB("coulombForce   default K = 0.00154629", coulombForce().GetK()==0.00154629);
  RowB("universalForce default G = 1.2331e-39", universalForce().GetG()==1.2331e-39);
  {
    coulombForce c2; c2.SetK(2*c2.GetK());
    pa.applyX(gX2 + 10.*gDir);
    Double_t rel = (c2.Evaluate(&pa,&pb) - 2*coul.Evaluate(&pa,&pb)).Mag() / coul.Evaluate(&pa,&pb).Mag();
    Row("coulombForce   SetK(2K) doubles the force", rel, 1e-15);
    universalForce g2; g2.SetG(3*g2.GetG());
    rel = (g2.Evaluate(&pa,&pb) - 3*grav.Evaluate(&pa,&pb)).Mag() / grav.Evaluate(&pa,&pb).Mag();
    Row("universalForce SetG(3G) triples the force", rel, 1e-15);
  }
}

// Alpha on Au, both moving, for trajectory comparisons.
void BuildPair(event & ev, force* f, Bool_t record){
  ev.AddForce(f);
  ev.AddParticle(new EMparticle(197, 79, TVector3(0.,0.,0.), TVector3(0.,0.,0.), false, record));
  ev.AddParticle(new EMparticle(4, 2, TVector3(-300., 20., 0.), TVector3(0.05, 0., 0.), false, record));
}

Bool_t SameTrajectory(force* fa, force* fb, const char* integ, Int_t n){
  event ea; BuildPair(ea, fa, kFALSE); ea.SetIntegrator(integ);
  event eb; BuildPair(eb, fb, kFALSE); eb.SetIntegrator(integ);
  ea.DeriveDTN(n, 1.);
  eb.DeriveDTN(n, 1.);
  Bool_t ok = kTRUE;
  Int_t i;
  for(i=0; i<2; i++){
    ok = ok && SameBits(ea.getParticle(i)->GetX(), eb.getParticle(i)->GetX());
    ok = ok && SameBits(ea.getParticle(i)->GetV(), eb.getParticle(i)->GetV());
  }
  return ok && SameBits(ea.GetNetE(), eb.GetNetE());
}

void CheckSoftening(){
  Header("5. softening");
  coulombForce   c0, cz;  cz.SetSoftening(0.);
  universalForce g0, gz;  gz.SetSoftening(0.);
  EMparticle pa(4, 2, TVector3(), TVector3(), false, false);
  EMparticle pb(197, 79, gX2, TVector3(), false, false);

  Double_t rs[] = {1e-3, 0.5, 7., 300., 1e5};
  Bool_t okC = kTRUE, okG = kTRUE;
  Int_t j;
  for(j=0; j<5; j++){
    pa.applyX(gX2 + rs[j]*gDir);
    okC = okC && SameBits(c0.Evaluate(&pa,&pb), cz.Evaluate(&pa,&pb)) && SameBits(c0.Potential(&pa,&pb), cz.Potential(&pa,&pb));
    okG = okG && SameBits(g0.Evaluate(&pa,&pb), gz.Evaluate(&pa,&pb)) && SameBits(g0.Potential(&pa,&pb), gz.Potential(&pa,&pb));
  }
  RowB("coulomb eps=0 : force, potential bit-identical to default", okC);
  RowB("gravity eps=0 : force, potential bit-identical to default", okG);
  RowB("coulomb eps=0 : 2000-step TAYLOR trajectory bit-identical", SameTrajectory(&c0, &cz, "TAYLOR", 2000));
  RowB("coulomb eps=0 : 2000-step VERLET trajectory bit-identical", SameTrajectory(&c0, &cz, "VERLET", 2000));
  RowB("coulomb eps=0 : 2000-step RK4 trajectory bit-identical",    SameTrajectory(&c0, &cz, "RK4", 2000));
  RowB("gravity eps=0 : 2000-step VERLET trajectory bit-identical", SameTrajectory(&g0, &gz, "VERLET", 2000));

  // eps > 0 : finite and bounded down to r = 0.
  // max over r of r / (r^2+eps^2)^(3/2) is 2/(3 sqrt(3) eps^2) at r = eps/sqrt(2).
  Double_t eps = 2.;
  coulombForce   cs; cs.SetSoftening(eps);
  universalForce gs; gs.SetSoftening(eps);
  Double_t kqq  = cs.GetK()*pa.GetQ()*pb.GetQ();
  Double_t gmm  = gs.GetG()*pa.GetM0()*pb.GetM0();
  Double_t fmax = 2./(3.*TMath::Sqrt(3.)*eps*eps);
  Double_t rz[] = {1e-2, 1e-6, 1e-12, 0.};
  Bool_t finC = kTRUE, finG = kTRUE;
  for(j=0; j<4; j++){
    pa.applyX(gX2 + rz[j]*gDir);
    Double_t fc = cs.Evaluate(&pa,&pb).Mag(), uc = cs.Potential(&pa,&pb);
    Double_t fg = gs.Evaluate(&pa,&pb).Mag(), ug = gs.Potential(&pa,&pb);
    finC = finC && TMath::Finite(fc) && TMath::Finite(uc) && fc<=kqq*fmax*(1+1e-12);
    finG = finG && TMath::Finite(fg) && TMath::Finite(ug) && fg<=gmm*fmax*(1+1e-12);
  }
  RowB("coulomb eps=2 : finite, |F| <= bound at r = 1e-2 ... 0", finC);
  RowB("gravity eps=2 : finite, |F| <= bound at r = 1e-2 ... 0", finG);
  pa.applyX(gX2);
  Row("coulomb eps=2 : |U(0) - kq1q2/eps| / U(0)", TMath::Abs(cs.Potential(&pa,&pb) - kqq/eps)/(kqq/eps), 1e-15);
  Row("gravity eps=2 : |U(0) + Gm1m2/eps| / U(0)", TMath::Abs(gs.Potential(&pa,&pb) + gmm/eps)/(gmm/eps), 1e-15);
  // far away the softened force goes to the plain one (relative difference ~ 1.5 eps^2/r^2)
  pa.applyX(gX2 + 1000.*gDir);
  Row("coulomb eps=2 at r=1000 : relative difference to eps=0",
      (cs.Evaluate(&pa,&pb) - c0.Evaluate(&pa,&pb)).Mag()/c0.Evaluate(&pa,&pb).Mag(), 1e-5);
}

void CheckRecordPeriod(){
  Header("6. record period (two moving recorded particles, dt = 1)");
  coulombForce coul;
  const Int_t N = 1000;
  const char* integs[] = {"TAYLOR", "VERLET"}; // applyForce() and finishStep() paths
  Int_t periods[] = {1, 7, 10, 1000, 0};       // 0 is taken as 1
  Int_t ii, jj;
  for(ii=0; ii<2; ii++){
    event ref; BuildPair(ref, &coul, kTRUE); ref.SetIntegrator(integs[ii]);
    ref.DeriveDTN(N, 1.);
    path* pr = ref.getParticle(1)->GetPath();
    for(jj=0; jj<5; jj++){
      Int_t n  = periods[jj];
      Int_t ne = (n<1) ? 1 : n;
      event ev; BuildPair(ev, &coul, kTRUE); ev.SetIntegrator(integs[ii]);
      ev.getParticle(0)->SetRecordPeriod(n);
      ev.getParticle(1)->SetRecordPeriod(n);
      ev.DeriveDTN(N, 1.);
      path* pp = ev.getParticle(1)->GetPath();
      Int_t expected = N/ne;
      // every kept point is the reference point of the same step, bit for bit
      Bool_t same = (pp->GetRecordSize()==expected) && (ev.getParticle(0)->GetPath()->GetRecordSize()==expected);
      Int_t k;
      for(k=0; same && k<pp->GetRecordSize(); k++){
        Int_t kr = (k+1)*ne - 1;
        same = SameBits(pp->GetRecordX(k), pr->GetRecordX(kr)) && SameBits(pp->GetRecordV(k), pr->GetRecordV(kr))
            && SameBits(pp->GetRecordT(k), pr->GetRecordT(kr));
      }
      // the trajectory itself does not depend on the record period
      same = same && SameBits(ev.getParticle(1)->GetX(), ref.getParticle(1)->GetX());
      RowB(TString::Format("%-6s n = %4d : %4d points of step n, 2n, ...", integs[ii], n, expected), same,
           TString::Format("%d points", pp->GetRecordSize()));
    }
  }
  { // an explicit AddRecord() is kept whatever the period
    EMparticle p(4, 2, TVector3(), TVector3(0.01,0,0), false, kTRUE);
    p.SetRecordPeriod(10);
    p.AddRecord();
    Int_t i;
    for(i=0; i<25; i++) p.applyForce(TVector3(), 1.);
    RowB("AddRecord() + 25 steps, n = 10 : 1 + 2 points", p.GetPath()->GetRecordSize()==3,
         TString::Format("%d points", p.GetPath()->GetRecordSize()));
  }
}

void CheckCopy(){
  Header("7. particle copy");
  coulombForce coul;
  event ev; BuildPair(ev, &coul, kTRUE); ev.SetIntegrator("VERLET");
  ev.getParticle(1)->SetRecordPeriod(3);
  ev.DeriveDTN(100, 1.);
  particle* src = ev.getParticle(1);
  particle cp(*src);
  RowB("copy keeps the step count", cp.GetECNT()==src->GetECNT() && cp.GetECNT()==100,
       TString::Format("%ld / %ld", cp.GetECNT(), src->GetECNT()));
  RowB("copy keeps the charge (GetQ on particle)", cp.GetQ()==2. && src->GetQ()==2.,
       TString::Format("%g", cp.GetQ()));
  RowB("copy keeps the record period and the path", cp.GetRecordPeriod()==3 &&
       cp.GetPath()->GetRecordSize()==src->GetPath()->GetRecordSize(),
       TString::Format("n = %d, %d points", cp.GetRecordPeriod(), cp.GetPath()->GetRecordSize()));
  RowB("copy keeps the time and the state", SameBits(cp.GetTime(), src->GetTime()) &&
       SameBits(cp.GetX(), src->GetX()) && SameBits(cp.GetV(), src->GetV()));
  // a copied event continues from the same step count
  event ec(ev);
  RowB("copied event keeps the step count of its particles", ec.getParticle(1)->GetECNT()==100 &&
       ec.getParticle(0)->GetECNT()==100 && ec.getParticle(1)->GetQ()==2. && ec.getParticle(0)->GetQ()==79.);
  // the copy of a derived event and the source continue identically
  ec.DeriveDTN(50, 1.);
  ev.DeriveDTN(50, 1.);
  RowB("copied event continues bit-identically to the source", SameBits(ec.getParticle(1)->GetX(), ev.getParticle(1)->GetX())
       && ec.getParticle(1)->GetPath()->GetRecordSize()==ev.getParticle(1)->GetPath()->GetRecordSize(),
       TString::Format("%d points", ec.getParticle(1)->GetPath()->GetRecordSize()));
}

void CheckHelpers(){
  Header("8. GetGamma, GetMomentum");
  particle p(4, TVector3(), TVector3(0.36, 0.48, 0.));   // |v| = 0.6 c, gamma = 1.25
  Row("|GetGamma() - 1.25|", TMath::Abs(p.GetGamma() - 1.25), 1e-15);
  RowB("GetGamma()*m0 == GetM(kFALSE)", SameBits(p.GetGamma()*p.GetM0(), p.GetM(kFALSE)));
  Row("|GetMomentum() - m0 v|", (p.GetMomentum() - 4.*p.GetV()).Mag(), 1e-15);
  RowB("GetMomentum() == GetP()", SameBits(p.GetMomentum(), p.GetP()));
  Row("|GetMomentum(kTRUE) - gamma m0 v|", (p.GetMomentum(kTRUE) - 1.25*4.*p.GetV()).Mag(), 1e-14);
}

int main(){
  CheckForces();
  CheckSoftening();
  CheckRecordPeriod();
  CheckCopy();
  CheckHelpers();

  std::cout<<std::endl<<(gFail==0 ? "ALL CHECKS PASSED" : TString::Format("%d CHECK(S) FAILED", gFail).Data())<<std::endl;
  return gFail==0 ? 0 : 1;
}
