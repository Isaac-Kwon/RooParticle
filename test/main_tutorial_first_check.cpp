////////////////////////////////////////////
//
// RooParticle tutorial, part 1 : the first event
//
// One alpha particle (m = 4 amu, q = 2 e) on a fixed gold-like target
// (q = 79 e), integrated for a fixed number of steps with velocity Verlet.
// Prints the final state, the scattering angle and the energy drift.
//
// The code between the //! [name] markers is shown in the manual
// (docs/tutorial.md). No arguments. Exit code is non-zero when a result is
// not what the tutorial says it is.
//
////////////////////////////////////////////

//! [includes]
#include "event.hpp"      // event, and through it particle, force, inspector
#include "particle.hpp"   // particle, EMparticle
#include "force.hpp"      // coulombForce and the other forces

#include "TMath.h"
#include "TVector3.h"
#include <cstdio>
//! [includes]

int main(){
  Int_t nfail = 0;

  //! [setup]
  // The force is not owned by the event: keep it alive as long as the event.
  coulombForce coulomb;

  event ev;
  ev.AddForce(&coulomb);

  // EMparticle(mass [amu], charge [e], position [fm], velocity [c], invincible, record)
  // The event owns the particles: create them with new and never delete them yourself.
  particle* target = ev.AddParticle(
      new EMparticle(197., 79., TVector3(0., 0., 0.), TVector3(0., 0., 0.), true, false));
  particle* alpha  = ev.AddParticle(
      new EMparticle(4., 2., TVector3(-1000., 50., 0.), TVector3(0.05, 0., 0.), false, false));
  //! [setup]

  //! [integrator]
  ev.SetIntegrator("VERLET");          // or ev.SetIntegrator(event::kVelocityVerlet)
  std::printf("integrator           : %s\n", ev.GetIntegratorName().Data());
  //! [integrator]

  //! [derive]
  // The alpha has to come in from 1000 fm and go out again to about 1000 fm:
  // about 2 x 1000 fm / 0.05 c = 40000 fm/c. With dt = 1 fm/c that is 40000 steps.
  const Double_t dt     = 1.;                                  // fm/c
  const Double_t T      = 2. * 1000. / alpha->GetV().Mag();    // fm/c
  const Int_t    nsteps = TMath::Nint(T / dt);

  const TVector3 vin = alpha->GetV();
  const Double_t E0  = ev.GetNetE();   // kinetic + potential energy [amu c^2]

  ev.DeriveDTN(nsteps, dt);            // nsteps calls of DeriveDT(dt)
  //! [derive]

  //! [result]
  TVector3 x = alpha->GetX();          // fm
  TVector3 v = alpha->GetV();          // c
  std::printf("steps, time          : %ld steps, %.1f fm/c\n", (long) alpha->GetECNT(), alpha->GetTime());
  std::printf("final position (fm)  : (%.3f, %.3f, %.3f)  r = %.3f fm\n", x.X(), x.Y(), x.Z(), x.Mag());
  std::printf("final velocity (c)   : (%.6f, %.6f, %.6f)  |v| = %.6f c\n", v.X(), v.Y(), v.Z(), v.Mag());

  // Scattering angle: between the initial and the final velocity.
  Double_t theta = vin.Angle(v) * TMath::RadToDeg();
  std::printf("scattering angle     : %.4f deg\n", theta);
  //! [result]

  //! [energy]
  const Double_t E1 = ev.GetNetE();
  const Double_t drift = TMath::Abs((E1 - E0) / E0);
  std::printf("energy  E0, E1       : %.12e  %.12e amu c^2\n", E0, E1);
  std::printf("relative drift       : %.2e\n", drift);
  std::printf("kinetic energy       : %.4f MeV at the start, %.4f MeV at the end\n",
              0.5 * 4. * vin.Mag2() * 931.494, alpha->GetKE() * 931.494);
  //! [energy]

  // ---- self-check (not shown in the manual)
  // The direction at r = 1000 fm is not yet the asymptotic one (see part 3),
  // so only a loose check on the angle; reference RK4, dt = 0.125 : 50.9 deg.
  if(!(TMath::Abs(theta - 50.9) < 0.3)){ std::printf("FAIL : angle\n"); nfail++; }
  if(!(drift < 1e-8))                   { std::printf("FAIL : energy drift\n"); nfail++; }
  if(alpha->GetECNT() != nsteps || TMath::Abs(alpha->GetTime() - T) > 1e-6){ std::printf("FAIL : steps / time\n"); nfail++; }
  if(target->GetX().Mag() != 0. || target->GetTime() != 0.){ std::printf("FAIL : target moved\n"); nfail++; }
  if(!(x.Mag() > 900. && x.Mag() < 1100.)){ std::printf("FAIL : final distance\n"); nfail++; }

  std::printf("%s\n", nfail == 0 ? "TUTORIAL 1 OK" : "TUTORIAL 1 FAILED");
  return nfail == 0 ? 0 : 1;
}
