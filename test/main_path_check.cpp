// Self-check of path::GetDCA (no simulation). Exit code != 0 on any failure.
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "path.hpp"

static int nFail = 0;

static void Check(const char* what, Bool_t ok){
  std::cout<<(ok ? "  [ OK ] " : "  [FAIL] ")<<what<<std::endl;
  if(!ok) nFail++;
}

static Bool_t Near(Double_t a, Double_t b, Double_t tol){
  return std::fabs(a-b)<tol;
}

// straight line x = x0 + v t sampled at t = 0, dt, 2dt, ...
static void Fill(path & p, TVector3 x0, TVector3 v, Int_t n, Double_t dt){
  for(Int_t i=0; i<n; i++){
    p.AddRecord(x0+v*(i*dt), v, i*dt);
  }
}

int main(){
  std::cout<<"path::GetDCA checks"<<std::endl;
  Double_t DCA, T;

  // 1) line along x at y=3 passing the origin: exact DCA 3 at t = 10.3
  //    (samples every 1.0, so the minimum lies between samples 10 and 11)
  std::cout<<"1) point DCA between two samples"<<std::endl;
  {
    path p;
    Fill(p, TVector3(-10.3,3.,0.), TVector3(1.,0.,0.), 21, 1.0);
    p.GetDCA(TVector3(0,0,0), DCA, T, false, false);
    Double_t dD = DCA, tD = T;
    std::cout<<"     discrete     DCA="<<dD<<" t="<<tD<<std::endl;
    p.GetDCA(TVector3(0,0,0), DCA, T);
    std::cout<<"     interpolated DCA="<<DCA<<" t="<<T<<" (exact 3, 10.3)"<<std::endl;
    Check("DCA is nonzero (argument is used)", DCA>1.);
    Check("interpolated DCA == 3", Near(DCA,3.,1e-9));
    Check("interpolated time == 10.3", Near(T,10.3,1e-9));
    Check("interpolated DCA closer to exact than discrete", std::fabs(DCA-3.)<std::fabs(dD-3.));
    Check("interpolated time closer to exact than discrete", std::fabs(T-10.3)<std::fabs(tD-10.3));
    Check("discrete time is a sample time", Near(tD,10.,1e-12));
  }

  // 2) two lines with known mutual DCA 4; separation along x is 20.8-2t,
  //    so the closest approach is at t = 10.4 (between samples)
  std::cout<<"2) path vs path DCA"<<std::endl;
  {
    path a, b;
    Fill(a, TVector3(0,0,0),    TVector3(1,0,0),  21, 1.0);
    Fill(b, TVector3(20.8,4,0), TVector3(-1,0,0), 21, 1.0);
    a.GetDCA(&b, DCA, T, false, false);
    Double_t dD = DCA, tD = T;
    a.GetDCA(&b, DCA, T);
    std::cout<<"     discrete DCA="<<dD<<" t="<<tD<<"  interpolated DCA="<<DCA<<" t="<<T<<" (exact 4, 10.4)"<<std::endl;
    Check("interpolated DCA == 4", Near(DCA,4.,1e-9));
    Check("interpolated time == 10.4", Near(T,10.4,1e-9));
    Check("interpolated better than discrete", std::fabs(DCA-4.)<std::fabs(dD-4.));

    // minimum at the path end falls back to the discrete value
    path c, d;
    Fill(c, TVector3(0,0,0), TVector3(1,0,0), 5, 1.0);
    Fill(d, TVector3(0,5,0), TVector3(0,0,0), 5, 1.0);
    c.GetDCA(&d, DCA, T);
    Check("minimum at path start: discrete value kept", Near(T,0.,1e-12) && Near(DCA,5.,1e-12));
  }

  // 3) failure cases report the sentinel (-1,-1) and do not crash
  std::cout<<"3) empty / mismatch -> sentinel (-1,-1), no crash"<<std::endl;
  {
    path e, p, q, r;
    Fill(p, TVector3(0,0,0), TVector3(1,0,0), 5, 1.0);
    Fill(q, TVector3(0,0,0), TVector3(1,0,0), 4, 1.0);
    Fill(r, TVector3(0,0,0), TVector3(1,0,0), 5, 0.5); // same length, other time grid
    DCA = 99.; T = 99.;
    e.GetDCA(TVector3(1,1,1), DCA, T);
    Check("empty path vs point", DCA==-1. && T==-1.);
    DCA = 99.; T = 99.;
    e.GetDCA(&p, DCA, T);
    Check("empty path vs non-empty path", DCA==-1. && T==-1.);
    DCA = 99.; T = 99.;
    p.GetDCA(&q, DCA, T);
    Check("length mismatch", DCA==-1. && T==-1.);
    DCA = 99.; T = 99.;
    e.GetDCA(&e, DCA, T);
    Check("two empty paths", DCA==-1. && T==-1.);
    DCA = 99.; T = 99.;
    p.GetDCA(&r, DCA, T);
    Check("time grid mismatch", DCA==-1. && T==-1.);
    DCA = 99.; T = 99.;
    p.GetDCA((path*)nullptr, DCA, T);
    Check("null pointer", DCA==-1. && T==-1.);
  }

  // 4) copy constructor accepts a const reference
  std::cout<<"4) copy constructor"<<std::endl;
  {
    path p;
    Fill(p, TVector3(0,0,0), TVector3(1,0,0), 3, 1.0);
    const path & cp = p;
    path q(cp);
    Check("copy has same size", q.GetMaxNumber()==3);
  }

  std::cout<<(nFail==0 ? "ALL PASSED" : "FAILURES")<<" ("<<nFail<<" failed)"<<std::endl;
  return nFail==0 ? 0 : 1;
}
