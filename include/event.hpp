#ifndef H_EVENT
#define H_EVENT 1;

#include "iostream"
#include "stdio.h"
#include "list"
#include "vector"

#include "TROOT.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"

#include "RPConfig.hpp"

#include "particle.hpp"
#include "force.hpp"
#include "eventVolume.hpp"
#include "inspector.hpp"
#include "recorder.hpp"

class inspector;

/// \brief One simulated system: particles, pairwise forces, inspectors and the time integration.
///
/// An event is filled with AddParticle(), AddForce() and AddInspector(), and is then
/// advanced with DeriveDT(), DeriveDTN(), DeriveMAX() or DeriveInspect(). Units are
/// fm, amu, c, e and time in fm/c (energies in amu c^2, momenta in amu c).
///
/// \code
/// event * ev = new event();
/// ev->AddForce(new coulombForce());
/// particle * target = ev->AddParticle(new EMparticle(206, 79, TVector3(0,0,0), TVector3(0,0,0), true, false));
/// particle * alpha  = ev->AddParticle(new EMparticle(4, 2, TVector3(-1000,50,0), TVector3(0.05,0,0), false, true));
/// ev->SetIntegrator("VERLET");
/// ev->AddInspector(new inspectorP(target, alpha, "DEG", 1));
/// ev->AddInspector(new inspectorP(target, alpha, "CNT", 50000));
/// ev->DeriveInspect();
/// \endcode
///
/// \par Forces
/// Every force added with AddForce() acts on every pair of particles of the event
/// (there is no way to restrict a force to some particles), and the forces are summed.
/// A force is evaluated as `Evaluate(p1, p2)`, the force on p1 by p2.
///
/// \par Many bodies
/// The simulator is built for three or more bodies, not only one projectile on a fixed
/// target. With N particles every pair interacts: the force on a particle is the sum over
/// all other particles and all forces, so the cost grows as N^2 per force evaluation. All
/// moving particles are moved together (except in oddevent). Inspectors on a pair
/// (inspectorP) look at two chosen particles only, while inspectors on the event
/// (inspectorE) look at all of them (net energies, smallest distance between any two
/// particles, largest acceleration, ...). The total momentum and angular momentum of the
/// moving particles are conserved only if no invincible particle is present and the
/// integrator treats action and reaction symmetrically.
///
/// \par Invincible particles
/// An invincible particle (see particle::IsInvincible()) acts on the others but is never
/// moved by any integrator and never receives a force. It is skipped in the kinetic
/// energy, momentum, angular momentum, acceleration and volume checks, but it is counted
/// in potential energies, in GetMinDistance() and in GetDeriveN(). It counts every step
/// (particle::GetECNT()) but its own time does not advance.
///
/// \par Ownership
/// The event owns, and its destructor deletes,
///  - the particles given to AddParticle(),
///  - the inspectors given to AddInspector(),
///  - the pre-derived copy made by preDerive().
///
/// So particles and inspectors must be made with `new`, and one of them must not be given
/// to two events (a copied event makes its own copies). The event does NOT own the forces
/// (AddForce()) nor the volume (SetVolume()): they are shared with the copies of the event
/// and may live on the stack. The pointers returned by AddParticle(), getParticle() and
/// getInspector() stay valid until the event is deleted.
///
/// \par Copying
/// `event(event &e)` makes a new event from a template event, to run many events with the
/// same setup. It copies the particles (state, step count, time, charge, path, as plain
/// particle objects), makes new inspectors of the same method, criterion and abort flag
/// that watch the new particles (inspectorP) or the new event (inspectorE), and copies the
/// settings: integrator, inspect mode, step cap, pair symmetry, relativistic flag, RK45
/// tolerances, and the volume and force pointers (shared, not duplicated). It does not copy
/// the RK45 suggested substep and statistics, the last stop reason, the pre-derived copy,
/// the Velocity-Verlet cache, the reference energy of the "EDG"/"EDL" inspectors, nor
/// inspectors that are neither inspectorP nor inspectorE (they are dropped). The copy of an
/// oddevent is a plain event with the default integrator, not an oddevent. The default
/// assignment is deleted.
///
/// \par Integrators and dynamics settings
/// The integrator is selected by SetIntegrator(), see the table at integrator. The default is
/// kTaylor. DeriveTaylor() (and so oddevent) ignores SetIntegrator(), SetRelativistic()
/// and SetPairSymmetric(); every other integrator honors them.
class event{
public:
  /// \brief Integration scheme used by DeriveDT().
  ///
  /// All of them evaluate every force before any particle is moved (oddevent does not).
  /// They are non-relativistic (rest mass) unless SetRelativistic(kTRUE) is called, which
  /// kTaylor ignores. The default is kTaylor, the original scheme.
  ///
  /// | Enumerator | Name | Scheme | Order |
  /// |---|---|---|---|
  /// | kTaylor | `"TAYLOR"` | `x += v dt + a dt^2/2` ; `v += a dt` (original scheme, default) | 1 |
  /// | kSymplecticEuler | `"SEULER"` | `v += a dt` ; `x += v dt` (symplectic) | 1 |
  /// | kVelocityVerlet | `"VERLET"` | `x += v dt + a dt^2/2` ; `v += (a + a') dt/2` (symplectic) | 2 |
  /// | kRK4 | `"RK4"` | classical Runge-Kutta | 4 |
  /// | kRK45 | `"RK45"` | Dormand-Prince 5(4), adaptive | 5 |
  ///
  /// kRK45: DeriveDT(dt) still advances exactly dt and counts one step; inside, it takes as
  /// many substeps as the tolerance needs (SetTolerance()) and remembers the last suggested
  /// substep for the next call. It costs 6 force evaluations per substep plus 1 per call, so
  /// give it a dt larger than the substep.
  enum integrator {
    kTaylor = 0,      ///< "TAYLOR", original scheme (default)
    kSymplecticEuler, ///< "SEULER", symplectic Euler
    kVelocityVerlet,  ///< "VERLET", velocity Verlet
    kRK4,             ///< "RK4", classical Runge-Kutta
    kRK45             ///< "RK45", adaptive Dormand-Prince
  };

  /// \brief How the non-abort inspectors are combined in Inspect().
  ///
  /// The AND default is the design: a typical condition is "the particle reached the
  /// asymptote AND went through a minimum number of steps".
  enum inspectMode {
    kInspectAND = 0, ///< stop when all of them are satisfied (default)
    kInspectOR       ///< stop when any of them is satisfied
  };

  /// \brief Why the last DeriveInspect() returned.
  enum stopReason {
    kStopNone = 0,   ///< DeriveInspect() was not called yet (or is running)
    kStopInspector,  ///< the inspectors were satisfied (normal end)
    kStopAbort,      ///< an abort inspector fired
    kStopMaxSteps,   ///< the safety cap of steps (SetMaxInspectSteps()) was reached
    kStopInvalid,    ///< an inspector has an undefined method; nothing was derived
    kStopFailed      ///< the integrator gave up (kRK45 on a singular force, see HasIntegratorFailed())
  };

  /// \brief Make an empty event: no particles, forces, inspectors nor volume.
  ///
  /// Defaults: kTaylor, AND mode, pair symmetry on, non-relativistic, step cap
  /// MAXINSPECTSTEPS, RK45 tolerances RK45RELTOL and RK45ABSTOL (see RPConfig.hpp).
  event();

  /// \brief Make a copy of an event, to run many events from one template.
  ///
  /// See "Copying" in the class description for what is copied (particles, inspectors,
  /// settings), what is shared (forces, volume) and what is not copied. The inspectors of
  /// the source must refer to particles of the source event; an inspectorP whose particle is
  /// not in the source event is linked to a wrong particle of the copy. The argument is a
  /// non-const reference, so a const event cannot be copied.
  event(event &e);

  /// \brief Delete the particles, the inspectors and the pre-derived copy (not the forces nor the volume).
  virtual ~event();

  /// \brief Deleted: the default assignment would copy the pointers and delete them twice.
  event & operator=(const event &) = delete;

  //Setup Methods

  /// \brief Add a particle to the event.
  /// \param particle_ particle made with `new`; the event takes the ownership and deletes it.
  /// \return the same pointer, for convenience.
  ///
  /// The particles are indexed by the order of addition (getParticle()).
  particle* AddParticle(particle* particle_){particles.push_back(particle_); nparticle++; return particle_;}

  /// \brief Add a force that acts on every pair of particles of the event.
  /// \param force_ force object; NOT owned and NOT copied by the event (it must outlive the
  ///        event, and is shared with copies of the event).
  ///
  /// Several forces are summed. Invalidates the Velocity-Verlet cache.
  void AddForce(force * force_){forces.push_back(force_); nforce++; accCacheX.clear();}

  /// \brief Add an inspector (a terminating condition for DeriveInspect()).
  /// \param inspector_ inspector made with `new`; the event takes the ownership and deletes it.
  void AddInspector(inspector* inspector_){inspectors.push_back(inspector_);}

  /// \brief Set the volume used by IsInsideVolume() and the "VOUT"/"VOUTA" inspectors.
  /// \param volume_ volume object; NOT owned (may be a null pointer, which means no volume).
  void SetVolume(eventVolume * volume_){volume = volume_;}

  /// \brief Select the integrator by enumerator (see integrator).
  /// \param integ_ new integrator. Invalidates the Velocity-Verlet cache.
  ///
  /// The RK45 suggested substep and statistics are kept.
  void       SetIntegrator(integrator integ_){integ = integ_; accCacheX.clear();}

  /// \brief Select the integrator by name.
  /// \param name "TAYLOR", "SEULER", "VERLET", "RK4" or "RK45" (case is ignored).
  /// \return kTRUE if the name was known. Otherwise a message goes to stderr, kFALSE is
  ///         returned and the integrator is not changed.
  Bool_t     SetIntegrator(TString name);

  /// \brief The integrator in use (default kTaylor).
  integrator GetIntegrator(){return integ;}

  /// \brief Name of the integrator in use, as accepted by SetIntegrator(TString).
  TString    GetIntegratorName();

  /// \brief Set the tolerances of kRK45.
  /// \param rel relative tolerance (default RK45RELTOL = 1e-8 in RPConfig.hpp)
  /// \param abs absolute tolerance, in fm for positions and in c for velocities (default
  ///        RK45ABSTOL = 1e-10)
  ///
  /// A substep is accepted when, for every moving particle, the estimated error is about
  /// below `abs + rel * |value|` (RMS over the moving particles and over position and
  /// velocity parts), where the value is the position (fm) and the velocity (c). In the
  /// relativistic mode the momentum divided by the rest mass is used instead of the velocity.
  /// Invalid values (negative, or both zero) are refused with a message on stderr and the
  /// previous tolerances are kept.
  void     SetTolerance(Double_t rel, Double_t abs);

  /// \brief Relative tolerance of kRK45.
  Double_t GetRelTolerance(){return rk45Rel;}

  /// \brief Absolute tolerance of kRK45 (fm for positions, c for velocities).
  Double_t GetAbsTolerance(){return rk45Abs;}

  /// \brief Substep suggested for the next kRK45 DeriveDT, in fm/c (0 : none yet).
  Double_t GetRK45Substep(){return rk45H;}

  /// \brief Number of accepted kRK45 substeps so far.
  Long_t   GetRK45NSubsteps(){return rk45NAccept;}

  /// \brief Number of rejected kRK45 substeps so far.
  Long_t   GetRK45NRejected(){return rk45NReject;}

  /// \brief Use Newton's third law to evaluate each pair once (default on).
  /// \param sym_ kTRUE : each pair is evaluated once, `Evaluate(p1, p2)`, and the opposite
  ///        force is given to p2. kFALSE : every ordered pair is evaluated; use it for a
  ///        force whose `Evaluate(p2, p1)` is not `-Evaluate(p1, p2)`.
  ///
  /// With many bodies the force on particle i is the sum over all j of `Evaluate(i, j)`; with
  /// the symmetry on, the pair (i, j) with i < j is evaluated once and its opposite is added
  /// to j, which halves the evaluations and makes the total force of the moving particles
  /// vanish exactly. Invincible particles never get a force either way. Used by kSymplecticEuler,
  /// kVelocityVerlet, kRK4, kRK45 and GetMaxAcceleration(). Ignored by kTaylor and by
  /// oddevent, which always evaluate every ordered pair. Invalidates the Velocity-Verlet cache.
  void   SetPairSymmetric(Bool_t sym_=kTRUE){pairSymmetric = sym_; accCacheX.clear();}

  /// \brief Whether the pair symmetry is in use (default kTRUE).
  Bool_t IsPairSymmetric(){return pairSymmetric;}

  /// \brief Switch the dynamics to special relativity (default off).
  /// \param rel_ kTRUE to switch on.
  ///
  /// With c = 1:
  /// \f[ \vec p = \gamma m_0 \vec v, \qquad \frac{d\vec p}{dt} = \vec F, \qquad
  ///     \vec v = \frac{\vec p}{\sqrt{m_0^2 + p^2}} \f]
  /// The integrators kSymplecticEuler, kVelocityVerlet, kRK4 and kRK45 then integrate the
  /// momentum. kTaylor (and oddevent) ignore it and stay non-relativistic. The speeds must
  /// stay below 1 (c); otherwise the result is NaN. Retardation (finite speed of the fields)
  /// is not included. Invalidates the Velocity-Verlet cache.
  ///
  /// Only some quantities follow the setting: GetNetMomentum() and GetNetAngularMomentum()
  /// use \f$\gamma m_0 v\f$ and the energy inspector "EDG"/"EDL" uses GetNetERel(). The
  /// others stay non-relativistic: GetNetKE() (default argument), GetNetE(), GetMaxAcceleration()
  /// (force / rest mass), the "KEG"/"KEL"/"PKG"/"PKL" inspectors, and particle::GetP().
  void   SetRelativistic(Bool_t rel_=kTRUE){relativistic = rel_; accCacheX.clear();}

  /// \brief Whether the relativistic flag is set (the flag only; kTaylor ignores it).
  Bool_t IsRelativistic(){return relativistic;}

  /// \brief Set how the non-abort inspectors are combined in Inspect() (default kInspectAND).
  void        SetInspectMode(inspectMode mode_){imode = mode_;}

  /// \brief The way the non-abort inspectors are combined.
  inspectMode GetInspectMode(){return imode;}

  /// \brief Set the safety cap of DeriveInspect().
  /// \param n maximum number of steps; 0 means unlimited. The default is MAXINSPECTSTEPS
  ///        (100,000,000, RPConfig.hpp).
  void        SetMaxInspectSteps(Long_t n){maxInspectSteps = n;}

  /// \brief The safety cap of steps of DeriveInspect() (0 : unlimited).
  Long_t      GetMaxInspectSteps(){return maxInspectSteps;}

  /// \brief Why the last DeriveInspect() ended (kStopNone before the first call).
  stopReason  GetStopReason(){return lastStop;}

  /// \brief Name of the last stop reason: "NONE", "INSPECTOR", "ABORT", "MAXSTEPS", "INVALID" or "FAILED".
  TString     GetStopReasonName();

  /// \brief True once the integrator has given up.
  ///
  /// kRK45 gives up when one call of DeriveDT() needs more than RK45MAXSUBSTEPS substeps,
  /// which happens when two bodies fall into an unsoftened attractive singularity. The
  /// particles are left at the last accepted state and DeriveInspect() ends with
  /// kStopFailed. The flag is never cleared: the state is not usable any more.
  Bool_t      HasIntegratorFailed(){return integFailed;}

  /// \brief Correct the speed of every non-invincible particle so that the given velocity is
  ///        the velocity far from the others.
  ///
  /// For every non-invincible particle the speed is changed along its present direction to
  /// \f$ v'^2 = v^2 - (U_{all} + U_{inv})/m_0 \f$, where \f$U_{all}\f$ is the potential energy of
  /// the particle with all the other particles and \f$U_{inv}\f$ the part with the invincible
  /// ones (see particle::calibrateVelocity()). So the interaction with an invincible particle
  /// is removed completely and the one with a moving particle by half (it is shared). Uses
  /// the rest mass and the non-relativistic kinetic energy. A particle at rest stays at rest.
  /// If the radicand is negative (the potential barrier is higher than the kinetic energy)
  /// the velocity becomes NaN.
  void CalibrateVelocity();

  //DerivingMethod

  /// \brief Advance the event by one time step with the selected integrator.
  /// \param dt time step in fm/c (a negative value runs backwards in time)
  ///
  /// Every non-invincible particle is moved and its time advances by dt. Every particle,
  /// invincible or not, counts one step (particle::GetECNT()) and records its path
  /// according to its record period. With kRK45 the call advances exactly dt and counts one
  /// step, but takes internal substeps. oddevent overrides this.
  virtual void DeriveDT(Double_t dt = 1);

  /// \brief Call DeriveDT() n times.
  /// \param n number of steps
  /// \param dt time step in fm/c
  void DeriveDTN(Int_t n = 1, Double_t dt = 1);

  /// \brief Call DeriveDT() MAXSTEPS times (100,000, RPConfig.hpp), ignoring the inspectors.
  /// \param dt time step in fm/c
  void DeriveMAX(Double_t dt = 1);

  /// \brief Derive until Inspect() is satisfied.
  /// \param iperiod the inspectors are checked every iperiod steps, the first time before
  ///        any step (values below 1 are taken as 1)
  /// \param dt time step in fm/c
  /// \return kTRUE on a normal end (the inspectors are satisfied), kFALSE otherwise; the
  ///         reason is given by GetStopReason().
  ///
  /// - The non-abort inspectors are combined by the inspect mode, AND by default
  ///   (SetInspectMode()). An abort inspector (inspector::SetAbort()) ends the event by
  ///   itself, whatever the others say (stop reason kStopAbort, return kFALSE).
  /// - An event without any inspector ends at once without a step (kStopInspector, kTRUE).
  ///   An event with only abort inspectors goes on until one of them fires.
  /// - The safety cap SetMaxInspectSteps() ends the event with kStopMaxSteps (return kFALSE)
  ///   after that many steps, if the inspectors are not satisfied before. The cap is counted
  ///   in steps, not in inspection periods.
  /// - If any inspector (abort or not) has an undefined method, nothing is derived, a
  ///   message goes to stderr and kStopInvalid is returned with kFALSE.
  /// - Because inspection is every iperiod steps, the event may overshoot the condition by
  ///   up to iperiod - 1 steps; abort inspectors are also checked only then.
  Bool_t DeriveInspect(Int_t iperiod=1, Double_t dt = 1);

  //GetMethod

  /// \brief Get a particle by the order of addition.
  /// \param index from 0 to getNParticle() - 1
  /// \return the particle (owned by the event), or a null pointer with a message on stderr
  ///         if index is too large. A negative index is not checked (it returns the last
  ///         particle, or is undefined for an empty event).
  particle * getParticle(Int_t index);

  /// \brief Get a force by the order of addition (index as in getParticle(); the force is not owned).
  force * getForce(Int_t index);

  /// \brief Get an inspector by the order of addition (index as in getParticle(); owned by the event).
  inspector * getInspector(Int_t index);

  /// \brief Number of particles.
  Int_t getNParticle(){return nparticle;}

  /// \brief Number of forces.
  Int_t getNForce(){return nforce;}

  /// \brief Number of inspectors.
  Int_t getNInspector(){return inspectors.size();}

  /// \brief Number of derived steps (particle::GetECNT()) over all particles, invincible ones included.
  /// \param minima kFALSE : the largest count ; kTRUE : the smallest count.
  /// \return the count; for an event without particles the extreme value of Long_t
  ///         (the "DXG"/"DXL"/"DNG"/"DNL" inspectors use it).
  Long_t GetDeriveN(const Bool_t minima=kFALSE);

  //Calculation

  /// \brief Potential energy of one particle with the other particles (amu c^2).
  /// \param index index of the particle (as in getParticle(); it must be valid, otherwise
  ///        the call dereferences a null pointer)
  /// \param OnlyInvincibles kTRUE : count only the interaction with the invincible particles
  /// \return sum over the other particles and all forces of `Potential(p, other)`.
  ///
  /// The energy of a pair appears in the value of both particles, so the sum over all
  /// particles counts every pair twice (GetNetPE() corrects it).
  Double_t GetPE(Int_t index, Bool_t OnlyInvincibles=kFALSE);

  /// \brief Same as GetPE(Int_t, Bool_t) for a particle pointer.
  /// \param p particle (it needs not be in the event; every particle of the event except p itself counts)
  /// \param OnlyInvincibles kTRUE : count only the interaction with the invincible particles
  Double_t GetPE(particle* p, Bool_t OnlyInvincibles=kFALSE);

  /// \brief Kinetic energy of the non-invincible particles (amu c^2).
  /// \param SR kFALSE (default) : \f$\frac12 m_0 v^2\f$ . kTRUE : \f$\frac12 \gamma m_0 v^2\f$,
  ///        which is not the relativistic kinetic energy and is not conserved
  ///        (use GetNetKERel()).
  /// \return the sum over the non-invincible particles.
  Double_t GetNetKE(const Bool_t SR=kFALSE);

  /// \brief Potential energy of the whole event (amu c^2).
  ///
  /// Every pair counts once, as half the sum of GetPE() over the particles (so the potential
  /// must be symmetric). The pairs with an invincible particle and between invincible
  /// particles are included (the latter is a constant).
  Double_t GetNetPE();

  /// \brief GetNetKE(SR) + GetNetPE(); conserved in the non-relativistic mode only (see GetNetERel()).
  Double_t GetNetE(const Bool_t SR=kFALSE){return GetNetKE(SR)+GetNetPE();}

  /// \brief Relativistic kinetic energy, the sum of \f$(\gamma-1) m_0\f$ over the non-invincible particles (amu c^2).
  ///
  /// GetNetKE(kTRUE) is \f$\frac12 \gamma m_0 v^2\f$, which is not conserved.
  Double_t GetNetKERel();

  /// \brief GetNetKERel() + GetNetPE(); the conserved energy of the relativistic mode.
  Double_t GetNetERel(){return GetNetKERel()+GetNetPE();}

  /// \brief Total momentum of the non-invincible particles (amu c).
  ///
  /// The momentum is the one of the integrated dynamics: \f$\gamma m_0 v\f$ when the
  /// relativistic mode is really in effect (SetRelativistic() with an integrator other than
  /// kTaylor, never in oddevent), \f$m_0 v\f$ otherwise. With an invincible particle it is
  /// not conserved, because the invincible particle takes the recoil.
  TVector3 GetNetMomentum();

  /// \brief Total angular momentum \f$\sum (\vec x - \vec o) \times \vec p\f$ of the non-invincible
  ///        particles (amu c fm).
  /// \param origin the point \f$\vec o\f$ in fm (default the origin)
  ///
  /// The momentum is the same as in GetNetMomentum().
  TVector3 GetNetAngularMomentum(TVector3 origin = TVector3());

  /// \brief The largest evolved time among the particles (fm/c). Invincible particles do not
  ///        advance their time, so it is 0 if all are invincible.
  Double_t GetTime();

  /// \brief The smallest distance between two particles, invincible ones included (fm).
  ///        The largest double for fewer than two particles.
  Double_t GetMinDistance();

  /// \brief The largest acceleration (force / rest mass) of the non-invincible particles (amu c^2/fm/amu = c^2/fm).
  ///
  /// Evaluates the forces at the present state without changing it. It uses the rest mass
  /// also in the relativistic mode, and honors SetPairSymmetric(). 0 if nobody moves.
  Double_t GetMaxAcceleration();

  /// \brief Whether the non-invincible particles are inside the volume (SetVolume()).
  /// \param all kTRUE : all of them are inside ; kFALSE : at least one is inside.
  /// \return kTRUE without a volume, and also when there is no non-invincible particle.
  Bool_t   IsInsideVolume(Bool_t all=kTRUE);


  //Inspection Deriving and Pre-Deriving.

  /// \brief Make a copy of the event and derive it until its inspectors are satisfied (pre-derivation).
  /// \todo Unfinished feature (initial position set-up by pre-derivation, see TODO.md).
  ///       The copy is private and has no accessor, so the only visible effect is the time
  ///       spent. A previous pre-derived copy is deleted first. Do not rely on it.
  void preDerive();

  /// \brief Placeholder for the set-up after preDerive().
  /// \todo Unfinished feature: the body is empty.
  void preDeriveSetup();

  /// \brief Delete the pre-derived copy.
  /// \return kTRUE if there was one, kFALSE otherwise.
  Bool_t resetPreDerive();

  /// \brief Evaluate the terminating condition once, without deriving.
  /// \return kTRUE if the event should stop.
  ///
  /// The abort inspectors are checked first: if one is satisfied it returns kTRUE at once
  /// and sets the stop reason to kStopAbort. Otherwise the non-abort inspectors are combined
  /// by the inspect mode (AND or OR). With only abort inspectors it returns kFALSE (go on);
  /// with no inspector at all it returns kTRUE. The stop reason is only ever set to
  /// kStopAbort here; it is reset by DeriveInspect().
  Bool_t Inspect();

  /// \brief Describe the event (particles and inspectors) as text.
  /// \param onlymechanic passed to particle::Print() (which ignores it)
  /// \param mute kFALSE (default) also prints the text to stdout
  /// \param pprint kTRUE puts line breaks
  /// \return the text.
  TString Print(Bool_t onlymechanic=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);

  /// \brief Set the initial condition from the present one.
  /// \warning Declared but NOT defined anywhere: calling it fails at link time.
  /// \todo Unfinished feature (see TODO.md).
  Bool_t SetInitial();

protected:
  std::vector<particle*> particles; ///< the particles, in the order of addition (owned)
  std::vector<force*> forces;       ///< the forces (not owned)

  Int_t nparticle = 0; ///< number of particles
  Int_t nforce = 0;    ///< number of forces

  //Integrators

  /// \brief The original scheme: all forces are evaluated for the present state of every
  ///        moving particle (every ordered pair), then every particle is moved.
  /// \param dt time step in fm/c
  ///
  /// Non-relativistic; ignores the pair symmetry and the relativistic flag.
  void DeriveTaylor(Double_t dt);

  /// \brief Symplectic Euler step, `v += a dt` then `x += v dt` (momentum in the relativistic mode).
  /// \param dt time step in fm/c
  void DeriveSymplecticEuler(Double_t dt);

  /// \brief Velocity-Verlet step (kick-drift-kick on the momentum in the relativistic mode).
  /// \param dt time step in fm/c
  ///
  /// The acceleration (the force in the relativistic mode) at the end of the step is reused
  /// at the start of the next one, as long as every particle is at the position where
  /// the last step left it. Only the positions are compared: changing a force parameter or a
  /// particle velocity between two steps is not noticed. AddForce(), SetIntegrator(),
  /// SetPairSymmetric() and SetRelativistic() clear the cache.
  void DeriveVelocityVerlet(Double_t dt);

  /// \brief Classical 4th order Runge-Kutta step.
  /// \param dt time step in fm/c
  void DeriveRK4(Double_t dt);

  /// \brief Dormand-Prince 5(4) with adaptive substeps, advancing exactly dt (see integrator).
  /// \param dt time step in fm/c ; may be negative
  ///
  /// The 5th order solution is kept and the embedded 4th order one gives the error. The
  /// substep suggested at the end of a call is used as the first substep of the next call.
  /// If the substep reaches 1e-12 of dt without meeting the tolerance (for example at a
  /// singular force), the step is accepted anyway with a warning on stderr (once per event).
  /// Counts one step for every particle per call, whatever the number of substeps.
  void DeriveRK45(Double_t dt);

  /// \brief Net force on every particle at the present state, in amu c^2/fm (zero for invincibles).
  /// \param frc output, resized to the number of particles
  void ComputeForce(std::vector<TVector3> & frc);

  /// \brief Acceleration F/m0 of every particle at the present state, in c^2/fm (zero for invincibles).
  /// \param acc output, resized to the number of particles
  void ComputeAcceleration(std::vector<TVector3> & acc);

  /// \brief Close one step: every particle counts the step, the moving ones advance their time, paths are recorded.
  /// \param dt time step in fm/c
  void FinishStep(Double_t dt);

  // The Runge-Kutta integrators work on a state (x, u) of every particle, with
  // u = v (non-relativistic) or u = p = gamma m0 v (relativistic).
  // Invincible particles are kept in the vectors but never changed.

  /// \brief Whether the relativistic option is really in effect: set, and not kTaylor.
  ///
  /// Virtual: oddevent overrides it with kFALSE because it never uses the option.
  virtual Bool_t UseRelativistic(){return relativistic && integ!=kTaylor;}

  /// \brief The state variable u of the velocity v: v itself, or \f$\gamma m_0 v\f$ when relativistic.
  TVector3 MomentumOf(particle * p, const TVector3 & v);

  /// \brief The velocity v of the state variable u (inverse of MomentumOf()).
  TVector3 VelocityOf(particle * p, const TVector3 & u);

  /// \brief Read the positions x and state variables u of all particles.
  void     GetState(std::vector<TVector3> & x, std::vector<TVector3> & u);

  /// \brief Write positions and state variables back to the non-invincible particles.
  void     SetState(const std::vector<TVector3> & x, const std::vector<TVector3> & u);

  /// \brief Time derivatives dx/dt = v and du/dt (F/m0, or F when relativistic) at the present state.
  void     ComputeDerivative(std::vector<TVector3> & dx, std::vector<TVector3> & du);

private:
  eventVolume * volume = nullptr; ///< not owned
  event * pEvent = nullptr;       ///< pre-derived copy (owned)
  Bool_t preDerived = kFALSE;     ///< whether pEvent is valid

  std::vector<inspector*> inspectors; ///< owned

  integrator  integ = kTaylor;            ///< selected integrator
  inspectMode imode = kInspectAND;        ///< how the inspectors are combined
  Long_t      maxInspectSteps = MAXINSPECTSTEPS; ///< safety cap of DeriveInspect()
  stopReason  lastStop = kStopNone;       ///< reason of the last DeriveInspect()

  Bool_t      pairSymmetric = kTRUE;  ///< Newton's third law in use
  Bool_t      relativistic  = kFALSE; ///< relativistic flag

  // kRK45 : tolerances, substep suggested for the next call, statistics
  Double_t    rk45Rel = RK45RELTOL; ///< relative tolerance
  Double_t    rk45Abs = RK45ABSTOL; ///< absolute tolerance
  Double_t    rk45H   = 0.;         ///< substep suggested for the next call
  Long_t      rk45NAccept = 0;      ///< accepted substeps
  Long_t      rk45NReject = 0;      ///< rejected substeps
  Bool_t      rk45Warned  = kFALSE; ///< the minimum-substep warning was printed
  Bool_t      integFailed = kFALSE; ///< the integrator gave up (see HasIntegratorFailed())

  // Velocity Verlet reuses the acceleration (the force in the relativistic mode)
  // of the previous step, as long as no particle was moved from outside in the meantime.
  std::vector<TVector3> accCache;  ///< acceleration (or force) at the end of the last step
  std::vector<TVector3> accCacheX; ///< positions for which accCache is valid
};

#endif
