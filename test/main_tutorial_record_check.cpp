////////////////////////////////////////////
//
// RooParticle tutorial, part 8 : recording tracks sparsely
//
// The alpha of part 3 recorded every step, every 10th and every 100th step.
// Shows the number of path points, what path::GetLast*() returns when the
// step count is not a multiple of the period, the initial point, and the
// distance of closest approach from a sparse path.
//
// No arguments. Exit code is non-zero when a check fails.
//
////////////////////////////////////////////

#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"
#include "inspector.hpp"
#include "path.hpp"
#include "asymptote.hpp"

#include "TMath.h"
#include "TVector3.h"
#include <cstdio>

int main(){
  Int_t nfail = 0;
  coulombForce coulomb;
  const Double_t kappa = coulomb.GetK() * 2. * 79.;
  std::printf("%7s %7s %8s %14s %14s %16s %12s\n", "period", "steps", "points", "GetLastT", "GetTime",
              "|LastX - X| fm", "DCA err (fm)");
  const Int_t periods[] = {1, 10, 100};
  for(Int_t k = 0; k < 3; k++){
    //! [record]
    event ev;
    ev.AddForce(&coulomb);
    ev.SetIntegrator("RK4");
    particle* t = ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
    particle* a = ev.AddParticle(new EMparticle(4., 2., TVector3(-1000., 50., 0.), TVector3(0.05, 0., 0.), false, true));
    a->SetRecordPeriod(periods[k]);  // keep the state after every n-th step
    a->AddRecord();                   // the initial state is not recorded otherwise
    ev.AddInspector(new inspectorP(t, a, "VRG", 0.));
    ev.AddInspector(new inspectorP(t, a, "DSTG", 1000.));
    ev.DeriveInspect();

    path* trk = a->GetPath();
    // GetLast*() is the last RECORDED step; the final state is GetX(), GetV(), GetTime().
    Double_t lag = (trk->GetLastX() - a->GetX()).Mag();
    Double_t dca, tdca;
    trk->GetDCA(t->GetX(), dca, tdca);   // parabola through the 3 samples around the minimum
    Double_t err = TMath::Abs(dca - asymptote(a, t, kappa).GetClosestApproach());
    std::printf("%7d %7ld %8d %14.1f %14.1f %16.4f %12.2e\n", periods[k], (long) a->GetECNT(),
                trk->GetRecordSize(), trk->GetLastT(), a->GetTime(), lag, err);
    //! [record]
    Int_t expect = 1 + (Int_t) (a->GetECNT() / periods[k]);
    if(trk->GetRecordSize() != expect) nfail++;
    if(trk->GetInitT() != 0. || trk->GetInitX() != TVector3(-1000., 50., 0.)) nfail++;
    if(periods[k] == 1 && lag != 0.) nfail++;
    if(periods[k] > 1 && a->GetECNT() % periods[k] != 0 && !(lag > 0.)) nfail++;
    if(!(err < (periods[k] == 1 ? 1e-6 : 1e-2))) nfail++;
  }
  std::printf("%s\n", nfail == 0 ? "TUTORIAL 8 OK" : "TUTORIAL 8 FAILED");
  return nfail == 0 ? 0 : 1;
}
