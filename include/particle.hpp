#ifndef H_PARTICLE
#define H_PARTICLE 1;

#include "TROOT.h"
#include "TString.h"
#include "TVector3.h"

// #include "event.hpp"
#include "path.hpp"

class inspector;
class event;

/// \brief A point body with mass, position, velocity, step count, time and an optional recorded path.
///
/// Units: mass in amu, position in fm, velocity in c, time in fm/c, momentum in amu c,
/// energy in amu c^2. The charge is 0 for a plain particle (see EMparticle).
///
/// \par Invincible particles
/// An invincible particle is a fixed one: it exerts forces on the others but is never moved
/// by the event (the integrators skip it and `applyForce()` leaves it alone). It still counts
/// every step (GetECNT(), so "CNT" inspectors work on it) and records its path, but its time
/// GetTime() does not advance. The `apply*` setters change it anyway when called directly.
///
/// \par Step count, time and path
/// Every derived step increments GetECNT(). A moving particle also advances its time by dt.
/// If the particle is recording, a point (x, v, time) is appended to the path after a step
/// when the step count is a multiple of the record period (SetRecordPeriod(), default 1).
/// The constructor does not record the initial state: call AddRecord() to store it.
///
/// \par Ownership and copying
/// The particle owns its path (deleted by the destructor, GetPath() only lends it). The copy
/// constructor copies the state, the step count and time, the charge, the record flag and
/// period, and the path (a deep copy). The default assignment operator is not suitable
/// (it would share the path pointer): use the copy constructor.
class particle{
public:
  // particle();

  /// \brief Make a particle.
  /// \param m0_ rest mass in amu
  /// \param x0_ initial position in fm
  /// \param v0_ initial velocity in c
  /// \param Invincible_ true for a fixed particle that acts on others but never moves
  /// \param Record_ true to keep the track of the particle in a path (GetPath()); the path is
  ///        then empty until the first recorded step
  ///
  /// The charge is 0 (use EMparticle for a charged particle).
  particle(Double_t m0_ = 1,
           //Double_t q_  = 1,
           TVector3 x0_ = TVector3(),
           TVector3 v0_ = TVector3(),
           Bool_t   Invincible_ = false,
           Bool_t   Record_     = false);

  //Force Methods (Holding and Release)

  /// \brief Copy constructor: state, step count, time, charge, record flag and period, and a copy of the path.
  /// \param p particle to copy
  ///
  /// The copy continues the same step count, so "CNT" inspectors count from there. A copy of
  /// an EMparticle is a plain particle with the same charge.
  particle(const particle & p);

  /// \brief Destructor; deletes the path (not virtual: delete a particle through its own or its base type only, no member is allocated by subclasses).
  ~particle(){if(Record)delete ppath;}

  /// \brief Add a force to the pending force of this particle.
  /// \param f force in amu c^2/fm. Several calls sum up until releaseForce() or resetHold().
  void holdForce(const TVector3 f);

  /// \brief Apply the held force for one step of dt (see applyForce()) and clear it.
  /// \param dt time step in fm/c
  ///
  /// This is the step of the original scheme: it counts a step also for an invincible
  /// particle (which then does not move) and also when no force was held.
  void releaseForce(const Double_t dt);

  /// \brief Discard the held force.
  void resetHold();

  /// \brief Close one integration step whose x, v were already updated by the event
  ///        (step count, time and path record). Used by the selectable integrators.
  /// \param dt time step in fm/c; the time advances by dt unless the particle is invincible.
  void finishStep(const Double_t dt);

  /// \brief Velocity change which makes the given velocity the one far from the other bodies.
  /// \param ev_ event the particle is in (the potentials come from its forces)
  /// \param ignoreInvincible false (default) : an invincible particle returns the zero vector;
  ///        true : it is calibrated like the others
  /// \return the vector to add to v (along the present direction of v, changing the speed
  ///         to \f$v'^2 = v^2 - (U_{all} + U_{inv})/m_0\f$, rest mass and non-relativistic).
  ///         NaN if the radicand is negative. The particle is not changed.
  ///
  /// \f$U_{all}\f$ is event::GetPE() of this particle with every other body and \f$U_{inv}\f$
  /// the part with invincible bodies, so an invincible partner is removed fully and a moving
  /// partner by half (the energy of that pair is shared by the two bodies).
  TVector3 calculateCalibrateVelocity(event* ev_, Bool_t ignoreInvincible=kFALSE);

  /// \brief Apply calculateCalibrateVelocity() to this particle.
  /// \param ev_ event the particle is in
  /// \param ignoreInvincible see calculateCalibrateVelocity()
  /// \param verbose true prints the kinetic energy before and after
  /// \return the velocity change that was added.
  TVector3 calibrateVelocity(event* ev_, Bool_t ignoreInvincible=kFALSE, Bool_t verbose=kFALSE);

  //Applying Methods
  // The apply* setters below change the state directly: no step is counted, no time
  // advances, nothing is recorded, and an invincible particle is not protected.

  /// \brief Set position and velocity (fm, c).
  void applyMechanic(const TVector3 x_, const TVector3 v_){x = x_; v = v_;}

  /// \brief Set the position (fm).
  void applyX(const TVector3 x_){x=x_;}

  /// \brief Set the velocity (c).
  void applyV(const TVector3 v_){v=v_;}

  /// \brief Move the particle for one step under a constant force (the original Taylor step).
  /// \param f force in amu c^2/fm
  /// \param dt time step in fm/c (default 1)
  /// \param SR false (default) : mass is the rest mass; true : mass is \f$\gamma m_0\f$
  ///        (a crude correction, not full relativistic dynamics)
  ///
  /// With a = f/m: `x += v dt + a dt^2/2`, `v += a dt`, time += dt. Always counts a step and
  /// records the path (according to the record period); an invincible particle is neither
  /// moved nor advanced in time.
  void applyForce(const TVector3 f, const Double_t dt=1., Bool_t SR=kFALSE);

  /// \brief Add dx (fm) to the position.
  void applyDX(const TVector3 dx){x += dx;}

  /// \brief Add dv (c) to the velocity.
  void applyDV(const TVector3 dv){v += dv;}

  /// \brief Add dx (fm) to the position and dv (c) to the velocity.
  void applyDelta(const TVector3 dx, const TVector3 dv){x += dx; v += dv;}

  /// \brief Set a name for printing.
  void SetName(const TString name_){name=name_;}

  /// \brief Append the present x, v and time to the path now, whatever the record period.
  ///
  /// Does nothing if the particle is not recording. Use it to store the initial state.
  void AddRecord(){if(Record)ppath->AddRecord(x,v,etime);}

  /// \brief Set the record period of the path.
  /// \param n a step is appended when the step count is a multiple of n (default 1 : every
  ///        step). n < 1 is taken as 1.
  ///
  /// With n > 1, path::GetLast*() is the last recorded step, NOT the final state of the
  /// particle unless the step count is a multiple of n: use GetX(), GetV(), GetTime() for
  /// the final state.
  void  SetRecordPeriod(Int_t n){recordPeriod = (n<1) ? 1 : n;}

  /// \brief The record period (see SetRecordPeriod()).
  Int_t GetRecordPeriod(){return recordPeriod;}

  /// \brief Position in fm (a copy).
  TVector3 GetX(){return TVector3(x);}

  /// \brief Velocity in c (a copy).
  TVector3 GetV(){return TVector3(v);}

  /// \brief Non-relativistic momentum \f$m_0 v\f$ in amu c (same as GetMomentum(kFALSE)).
  ///
  /// Not the momentum of the relativistic dynamics, see GetMomentum().
  TVector3 GetP(){return TVector3(v * GetM());}

  /// \brief Momentum in amu c.
  /// \param relativistic false (default) : \f$m_0 v\f$ ; true : \f$\gamma m_0 v\f$
  TVector3 GetMomentum(Bool_t relativistic=kFALSE);

  /// \brief Number of derived steps ("evolved count"), counted also for an invincible particle.
  Long_t   GetECNT(){return ecount;}

  /// \brief Rest mass in amu.
  Double_t GetM0(){return m0;}

  /// \brief Mass in amu.
  /// \param restmass kTRUE (default) : the rest mass; kFALSE : \f$\gamma m_0\f$ (NaN for speed >= c).
  Double_t GetM(const Bool_t restmass=kTRUE);

  /// \brief Lorentz factor \f$1/\sqrt{1 - v^2}\f$ (v in c); NaN for speed above c.
  Double_t GetGamma();

  /// \brief Kinetic energy in amu c^2.
  /// \param SR kFALSE (default) : \f$\frac12 m_0 v^2\f$ ; kTRUE : \f$\frac12 \gamma m_0 v^2\f$.
  ///
  /// Neither is the relativistic kinetic energy \f$(\gamma-1) m_0\f$; for that see
  /// event::GetNetKERel(). Not zero for an invincible particle (the event skips it).
  Double_t GetKE(const Bool_t SR=kFALSE){return 0.5*GetM(!SR)*(v.Mag2());}

  /// \brief Charge in e. 0 for a plain particle, set by EMparticle.
  ///
  /// It lives here (not only in EMparticle) so that forces read it from particle* without
  /// casting, and so that it survives event copies, which copy particles as particle.
  Double_t GetQ(){return q;}

  /// \brief Evolved time in fm/c; does not advance for an invincible particle.
  Double_t GetTime(){return etime;}

  /// \brief The recorded path (owned by the particle; a null pointer if not recording).
  path*    GetPath(){return ppath;}

  /// \brief Whether the particle is fixed (invincible).
  Bool_t   IsInvincible(){return Invincible;}

  /// \brief Whether the particle records its path.
  Bool_t   IsRecording(){return Record;}

  /// \brief Describe the particle as text "[ Particle x: [..], v: [..] ]".
  /// \param onlymechanic has no effect at present
  /// \param mute false (default) also prints the text to stdout
  /// \return the text.
  TString  Print(Bool_t onlymechanic=kFALSE, Bool_t mute=kFALSE);

  friend class inspectorP;

protected:
  Double_t q          = 0; ///< charge in e

private:
  TString  name;                ///< name for printing
  TVector3 x          = TVector3(); ///< position in fm
  TVector3 v          = TVector3(); ///< velocity in c
  Double_t m0         = 1;      ///< rest mass in amu

  Bool_t   Invincible = false;  ///< fixed particle
  Bool_t   hold       = false;  ///< a force is held
  TVector3 hF         = TVector3(); ///< held force in amu c^2/fm

  Double_t etime      = 0.;     ///< evolved time in fm/c
  Long_t   ecount      = 0;     ///< number of derived steps

  Bool_t   Record     = false;  ///< records the path
  Int_t    recordPeriod = 1;    ///< a step is recorded when ecount is a multiple of it
  path*    ppath      = nullptr; ///< the path (owned)

  // Path record at the end of a step, following the record period.
  void recordStep(){if(Record && ecount%recordPeriod==0)ppath->AddRecord(x,v,etime);}
};

/// \brief A particle with a charge, for electromagnetic forces (coulombForce).
///
/// Charge in e, mass in amu. The charge is stored in particle::GetQ(). Every other
/// behavior is that of particle.
class EMparticle : public particle{
public:
  // EMparticle(): particle(){};

  /// \brief Make a charged particle.
  /// \param m0_ rest mass in amu
  /// \param q_ charge in e (default 1)
  /// \param x0_ initial position in fm
  /// \param v0_ initial velocity in c
  /// \param Invincible_ true for a fixed particle (for example the target nucleus)
  /// \param Record_ true to keep the path
  EMparticle(Double_t m0_ = 1,
             Double_t q_  = 1,
             TVector3 x0_ = TVector3(),
             TVector3 v0_ = TVector3(),
             Bool_t   Invincible_ = false,
             Bool_t   Record_     = false):
             particle(m0_, x0_, v0_, Invincible_, Record_){
               q = q_;
             }
  ~EMparticle(){}
  // GetQ() is inherited from particle.


private:
  // Double_t q;
};

#endif
