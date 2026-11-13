////////////////////////////////////////////
//
// RooParticle tutorial, part 3 : getting the physics out (two bodies)
//
// Alpha on a fixed gold-like target, b = 50 fm, v = 0.05 c, start at
// x = -1000 fm. After the event:
//   - the raw scattering angle from the final velocity,
//   - the exact asymptotic angle, speed and impact parameter (asymptote),
//     compared with the Rutherford formula at nominal and asymptotic values,
//   - the distance of closest approach from the recorded path.
// Also: event::CalibrateVelocity() to make the nominal speed the speed at infinity.
//
// No arguments. Exit code is non-zero when a number is off.
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

static Int_t gFail = 0;
void Expect(const char* what, Bool_t ok){
  if(!ok){ std::printf("  FAIL : %s\n", what); gFail++; }
}

int main(){
  const Double_t b0 = 50., v0 = 0.05;   // nominal impact parameter [fm] and speed [c]

  //! [run]
  coulombForce coulomb;
  event ev;
  ev.AddForce(&coulomb);
  particle* target = ev.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  // record = true : the alpha keeps its track (position, velocity, time of every step)
  particle* alpha  = ev.AddParticle(
      new EMparticle(4., 2., TVector3(-1000., b0, 0.), TVector3(v0, 0., 0.), false, true));
  ev.SetIntegrator("VERLET");
  ev.AddInspector(new inspectorP(target, alpha, "DEG", 1.));
  ev.AddInspector(new inspectorP(target, alpha, "CNT", 50000));

  // kappa of U = kappa / r : k q1 q2
  const Double_t kappa = coulomb.GetK() * alpha->GetQ() * target->GetQ();
  asymptote before(alpha, target, kappa);   // from the initial state, before any step

  ev.DeriveInspect();
  //! [run]

  //! [raw]
  // Raw angle: initial velocity vs velocity at the stop (r = 2821 fm here)
  Double_t thetaRaw = TVector3(v0, 0., 0.).Angle(alpha->GetV());
  std::printf("raw angle at r = %.0f fm          : %.6f deg\n",
              (alpha->GetX() - target->GetX()).Mag(), thetaRaw * TMath::RadToDeg());
  //! [raw]

  //! [asymptote]
  // asymptote(p1, p2, kappa): conserved energy, angular momentum and
  // Laplace-Runge-Lenz vector of the pair -> the state at t = -inf and +inf.
  asymptote after(alpha, target, kappa);
  std::printf("status                          : %s\n", after.GetStatusName().Data());
  std::printf("thetaCM, extrapolated at stop   : %.6f deg\n", after.GetThetaCM() * TMath::RadToDeg());
  std::printf("thetaCM, from the initial state : %.6f deg\n", before.GetThetaCM() * TMath::RadToDeg());
  std::printf("v_inf                           : %.7f c   (nominal %.2f c)\n", after.GetVinf(), v0);
  std::printf("effective b = |L| / (mu v_inf)  : %.4f fm  (nominal %.1f fm)\n",
              after.GetImpactParameter(), b0);
  std::printf("tilt of the incoming asymptote  : %.4f deg\n",
              asymptote::Angle(after.GetVInRel(), TVector3(1., 0., 0.)) * TMath::RadToDeg());
  std::printf("lab angle vs the x axis at +inf : %.6f deg\n",
              asymptote::Angle(TVector3(1., 0., 0.), after.GetV1Out()) * TMath::RadToDeg());

  const Double_t mu = after.GetMu();   // the alpha mass, the target is fixed
  std::printf("Rutherford at nominal (b, v)    : %.6f deg\n",
              asymptote::RutherfordAngle(b0, v0, mu, kappa) * TMath::RadToDeg());
  std::printf("Rutherford at (b_eff, v_inf)    : %.6f deg\n",
              asymptote::RutherfordAngle(after.GetImpactParameter(), after.GetVinf(), mu, kappa)
              * TMath::RadToDeg());
  //! [asymptote]

  //! [dca]
  // Closest approach to the (fixed) target position from the recorded path:
  // discrete minimum over the steps, then refined by a parabola through
  // the three samples around it (default).
  Double_t dca, tdca, dcaRaw, tRaw;
  alpha->GetPath()->GetDCA(target->GetX(), dcaRaw, tRaw, false, false);  // no interpolation
  alpha->GetPath()->GetDCA(target->GetX(), dca, tdca);                    // interpolated
  std::printf("DCA, discrete                   : %.9f fm at t = %.1f fm/c\n", dcaRaw, tRaw);
  std::printf("DCA, interpolated               : %.9f fm at t = %.3f fm/c\n", dca, tdca);
  std::printf("closest approach, asymptote     : %.9f fm\n", after.GetClosestApproach());
  std::printf("path points                     : %d\n", alpha->GetPath()->GetRecordSize());
  //! [dca]

  //! [calibrate]
  // CalibrateVelocity(): rescale the speed of every moving particle so that the
  // total energy equals the nominal kinetic energy (direction unchanged).
  // Then the nominal speed is the speed at infinity; b is still not b_eff.
  event cal;
  cal.AddForce(&coulomb);
  particle* t2 = cal.AddParticle(new EMparticle(197., 79., TVector3(), TVector3(), true, false));
  particle* a2 = cal.AddParticle(new EMparticle(4., 2., TVector3(-1000., b0, 0.), TVector3(v0, 0., 0.), false, false));
  cal.CalibrateVelocity();
  asymptote ac(a2, t2, kappa);
  std::printf("after CalibrateVelocity         : |v| = %.7f c, v_inf = %.9f c, b_eff = %.4f fm\n",
              a2->GetV().Mag(), ac.GetVinf(), ac.GetImpactParameter());
  //! [calibrate]

  // ---- self-check
  Double_t thAsy = after.GetThetaCM();
  Expect("status unbound", after.GetStatus() == asymptote::kUnbound);
  Expect("asymptote invariant along the run", TMath::Abs(thAsy - before.GetThetaCM()) < 1e-9);
  Expect("asymptote = Rutherford(b_eff, v_inf)",
         TMath::Abs(thAsy - asymptote::RutherfordAngle(after.GetImpactParameter(), after.GetVinf(), mu, kappa)) < 1e-12);
  Expect("nominal Rutherford is off by about 1 deg",
         TMath::Abs(asymptote::RutherfordAngle(b0, v0, mu, kappa) - thAsy) * TMath::RadToDeg() > 0.5);
  Expect("raw angle off by less than 0.1 deg", TMath::Abs(thetaRaw - thAsy) * TMath::RadToDeg() < 0.1);
  Expect("interpolated DCA within 1e-6 fm", TMath::Abs(dca - after.GetClosestApproach()) < 1e-6);
  // the parabola minimum lies between the samples, below the smallest one
  Expect("interpolated minimum below the discrete one", dca <= dcaRaw && TMath::Abs(tdca - tRaw) < 1.);
  Expect("calibrated v_inf = nominal", TMath::Abs(ac.GetVinf() - v0) < 1e-12);

  std::printf("%s\n", gFail == 0 ? "TUTORIAL 3 OK" : "TUTORIAL 3 FAILED");
  return gFail == 0 ? 0 : 1;
}
