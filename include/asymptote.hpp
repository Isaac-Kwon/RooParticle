/// \file asymptote.hpp
/// \brief Analytic two-body 1/r asymptote calculator (\ref asymptote), used to validate the integrators.

#ifndef H_ASYMPTOTE
#define H_ASYMPTOTE 1

#include "TROOT.h"
#include "TString.h"
#include "TVector3.h"

class particle;

/// \brief Analytic two-body 1/r asymptote calculator.
///
/// The simulation starts and stops at finite distance, where the potential
/// energy is not negligible. For two bodies interacting through a pure 1/r
/// potential the relative orbit is a conic section, and the state at infinity
/// follows exactly from the state at ANY single instant through the conserved
/// orbital energy, angular momentum and Laplace-Runge-Lenz vector.
/// asymptote takes that one state and returns the incoming and outgoing
/// asymptotes (speed at infinity, directions, scattering angle, effective
/// impact parameter, closest approach, lab velocities).
///
/// \warning TWO BODIES ONLY. RooParticle is built for three or more bodies,
/// but this class is exact only for exactly two interacting bodies (plus
/// nothing else): with a third body present the result is merely the
/// two-body orbit through the given instant and can be badly wrong. Do not use
/// it to predict the outcome of an N-body event. An N-body counterpart is
/// being added separately; until then this class serves as the analytic
/// reference for the two-body validation runs.
///
/// \section asy_pot Potential and units
/// \f[ U(r) = \frac{\kappa}{r} \f]
/// \f$\kappa > 0\f$ is repulsive, \f$\kappa < 0\f$ attractive. For the Coulomb
/// force \f$\kappa = k\,q_1 q_2\f$ with
/// \f$k = 0.00154629\ \mathrm{amu\,fm\,c^2/e^2}\f$. Units are fm, amu, c,
/// fm/c (any consistent set works: nothing is hard-coded). Angles are in
/// radians.
///
/// \section asy_valid Validity
/// The result is EXACT only when all of these hold:
///  - exactly two bodies (any third body, or a background field, breaks it),
///  - a pure 1/r interaction (no screening, no nuclear force, no softening),
///  - non-relativistic mechanics (rest masses, Newton's law; the relativistic
///    orbit is a precessing rosette, and the error grows like \f$v^2\f$),
///  - the given state is itself correct (an integration error already made
///    is propagated as is).
///
/// Otherwise it is only an approximation that improves as the given state is
/// taken farther from the other bodies.
///
/// \section asy_frames Frames
/// Relative coordinates are those of particle 1 with respect to particle 2:
/// \f$\vec r = \vec x_1 - \vec x_2\f$, \f$\vec v = \vec v_1 - \vec v_2\f$.
///  - Both free: reduced mass \f$\mu = m_1 m_2/(m_1+m_2)\f$, relative
///    coordinates, and the (constant) centre-of-mass velocity \f$\vec V_{cm}\f$.
///  - One fixed (an invincible particle, a fixed centre): \f$\mu\f$ is the
///    mass of the moving particle, the fixed one counts as at rest (its
///    velocity is ignored, as the event never moves it), lab coordinates are
///    the relative coordinates, \f$\vec V_{cm} = 0\f$.
///  - Both fixed: status \ref asymptote::kInvalid.
///
/// The masses taken from a \ref particle are its rest masses (particle::GetM0),
/// even if the event integrates relativistically.
///
/// \section asy_elements Orbital elements (of the relative motion)
/// \f[ E = \frac{\mu v^2}{2} + \frac{\kappa}{r} \f]
/// \f[ \vec L = \mu\,\vec r \times \vec v \f]
/// \f[ \vec A = (\mu \vec v) \times \vec L + \mu\kappa\,\hat r \f]
/// \f$E\f$ is the orbital energy, \f$\vec L\f$ the angular momentum and
/// \f$\vec A\f$ the Laplace-Runge-Lenz vector. For either sign of
/// \f$\kappa\f$, \f$\vec A\f$ points from the centre toward the point of
/// closest approach (periapsis), and \f$|\vec A| = \mu|\kappa|\,e\f$. The
/// eccentricity vector is defined here as
/// \f$\vec e = \vec A/(\mu|\kappa|)\f$, so it also points toward periapsis
/// (for \f$\kappa > 0\f$ this is minus the textbook \f$\vec e\f$ vector of a
/// repulsive orbit).
///
/// \section asy_status Status
/// See \ref asymptote::status. Only kUnbound, kHeadOn and kFree have
/// asymptotes (HasAsymptote()). For every other state the asymptotic getters
/// (speed at infinity, directions, angles, b, lab velocities) return 0 or zero
/// vectors, never NaN. What remains set for each state:
///  - kUnbound, kHeadOn, kFree: everything.
///  - kRadialAttractive, kParabolic, kBound: energy, L, A, eccentricity,
///    periapsis direction, closest approach; the asymptotic quantities and
///    lab velocities are zero.
///  - kInvalid: only what was set before the failure (the inputs kappa, masses
///    and flags; r, v and mu once they were computed); everything else is zero.
///    If the pointer-based constructor gets a null or identical pointer, even
///    the inputs are zero.
///
/// \section asy_example Example
/// \code
/// asymptote asy(alpha, target, 0.00154629*2*79);  // any instant of the event
/// if(asy.HasAsymptote()){
///   Double_t theta = asy.GetThetaCM()*TMath::RadToDeg();
///   TVector3 vout  = asy.GetV1Out();               // alpha at infinity (lab)
///   Double_t b     = asy.GetImpactParameter();     // effective b
/// }
/// \endcode
class asymptote{
public:
  /// \brief Classification of the relative orbit. Only kUnbound, kHeadOn and
  /// kFree have asymptotes (see HasAsymptote()).
  enum status {
    kUnbound = 0,    ///< \f$E>0\f$, \f$L\neq 0\f$: hyperbola, every asymptotic quantity is set
    kHeadOn,         ///< \f$E>0\f$, \f$L=0\f$, \f$\kappa>0\f$: straight back-scattering, \f$\theta=180^\circ\f$, \f$r_{min}=\kappa/E\f$, \f$b=0\f$
    kFree,           ///< \f$\kappa=0\f$: straight line, \f$\theta=0\f$, \f$r_{min}=b\f$ (with the relative velocity zero too it is reported as kParabolic with \f$r_{min}=r\f$)
    kRadialAttractive, ///< \f$L=0\f$, \f$\kappa<0\f$ (whatever the sign of E): falls into the singularity, no asymptote (\f$r_{min}=0\f$)
    kParabolic,      ///< \f$E=0\f$ (within round-off): \f$v_\infty=0\f$, no direction
    kBound,          ///< \f$E<0\f$: ellipse, no asymptote (E, L, A, e, \f$r_{min}\f$ are set)
    kInvalid         ///< zero distance, non-positive mass (of a moving particle), both particles fixed, or null/identical particle pointers
  };

  /// \brief Empty object; the status is kInvalid until Compute() is called.
  asymptote(){}
  /// \brief From two particles of an event (IsInvincible() selects the fixed centre).
  ///
  /// Uses the current position, velocity and rest mass of each particle.
  /// \param p1 particle 1 (the relative coordinate is \f$\vec x_1-\vec x_2\f$); not owned
  /// \param p2 particle 2; not owned
  /// \param kappa interaction constant of \f$U=\kappa/r\f$ (amu fm c^2 for the Coulomb force: k q1 q2)
  asymptote(particle* p1, particle* p2, Double_t kappa);
  /// \brief From raw states. \c fixed1 / \c fixed2: that particle is a fixed centre.
  /// \param x1 position of particle 1 (fm)
  /// \param v1 velocity of particle 1 (c); ignored if \c fixed1
  /// \param m1 mass of particle 1 (amu); only needs to be positive if not fixed
  /// \param fixed1 particle 1 is a fixed centre
  /// \param x2 position of particle 2 (fm)
  /// \param v2 velocity of particle 2 (c); ignored if \c fixed2
  /// \param m2 mass of particle 2 (amu); only needs to be positive if not fixed
  /// \param fixed2 particle 2 is a fixed centre
  /// \param kappa interaction constant
  asymptote(TVector3 x1, TVector3 v1, Double_t m1, Bool_t fixed1,
            TVector3 x2, TVector3 v2, Double_t m2, Bool_t fixed2,
            Double_t kappa);
  /// \brief From the relative state only.
  ///
  /// \f$\vec r=\vec x_1-\vec x_2\f$, \f$\vec v=\vec v_1-\vec v_2\f$, reduced
  /// mass \f$\mu\f$. Same as a moving particle 1 of mass \f$\mu\f$ around a
  /// fixed particle 2 at the origin (so IsFixedCentre() is true, the centre-of-mass
  /// velocity is zero and the lab velocities are the relative ones).
  /// \param r relative position (fm)
  /// \param v relative velocity (c)
  /// \param mu reduced mass (amu)
  /// \param kappa interaction constant
  asymptote(TVector3 r, TVector3 v, Double_t mu, Double_t kappa);

  /// \brief Recompute from two particles (arguments as the matching constructor); all previous results are discarded.
  /// \return HasAsymptote()
  Bool_t Compute(particle* p1, particle* p2, Double_t kappa);
  /// \brief Recompute from raw states (arguments as the matching constructor); all previous results are discarded.
  /// \return HasAsymptote()
  Bool_t Compute(TVector3 x1, TVector3 v1, Double_t m1, Bool_t fixed1,
                 TVector3 x2, TVector3 v2, Double_t m2, Bool_t fixed2,
                 Double_t kappa);
  /// \brief Recompute from the relative state (arguments as the matching constructor); all previous results are discarded.
  /// \return HasAsymptote()
  Bool_t Compute(TVector3 r, TVector3 v, Double_t mu, Double_t kappa);

  /// \name Status
  /// @{
  status   GetStatus()      const {return st;} ///< classification of the orbit
  /// \brief Name of the status as text: "UNBOUND", "HEADON", "FREE", "RADIAL_ATTRACTIVE", "PARABOLIC", "BOUND" or "INVALID".
  TString  GetStatusName()  const;
  Bool_t   HasAsymptote()   const {return st==kUnbound || st==kHeadOn || st==kFree;} ///< true for kUnbound, kHeadOn, kFree
  Bool_t   IsFixedCentre()  const {return fixedCentre;} ///< exactly one of the two bodies is fixed
  Bool_t   IsOutgoing()     const {return r.Dot(v) > 0.;} ///< the given state is past closest approach (\f$\vec r\cdot\vec v>0\f$); no check of the status
  /// @}

  /// \name Input, relative motion
  /// @{
  Double_t GetKappa()       const {return kappa;} ///< interaction constant of \f$U=\kappa/r\f$
  Double_t GetMu()          const {return mu;}         ///< reduced (or moving) mass (amu)
  TVector3 GetR()           const {return r;}          ///< \f$\vec x_1-\vec x_2\f$ (fm)
  TVector3 GetV()           const {return v;}          ///< \f$\vec v_1-\vec v_2\f$ (c)
  TVector3 GetVcm()         const {return vcm;}        ///< centre-of-mass velocity (c; 0 for a fixed centre)
  /// @}

  /// \name Conserved quantities and orbital elements
  /// @{
  Double_t GetEnergy()      const {return energy;}     ///< \f$E=\mu v^2/2+\kappa/r\f$ (amu c^2)
  TVector3 GetL()           const {return lvec;}       ///< \f$\mu\,\vec r\times\vec v\f$
  Double_t GetLMag()        const {return lvec.Mag();} ///< \f$|\vec L|\f$
  TVector3 GetLRL()         const {return avec;}       ///< Laplace-Runge-Lenz vector \f$(\mu\vec v)\times\vec L+\mu\kappa\,\hat r\f$
  TVector3 GetEccentricityVector() const {return evec;} ///< toward periapsis, \f$|\vec e|=e\f$
  Double_t GetEccentricity() const {return ecc;}       ///< eccentricity; 0 for kFree / kInvalid
  TVector3 GetPeriapsisDirection() const {return phat;}///< unit vector, centre to closest approach
  Double_t GetClosestApproach() const {return rmin;}   ///< distance of closest approach (fm)
  /// @}

  /// \name Asymptotes (relative motion)
  /// @{
  Double_t GetVinf()        const {return vinf;}       ///< relative speed at infinity (c)
  TVector3 GetVInRel()      const {return vinRel;}     ///< incoming relative velocity at \f$t\to-\infty\f$
  TVector3 GetVOutRel()     const {return voutRel;}    ///< outgoing relative velocity at \f$t\to+\infty\f$
  Double_t GetThetaCM()     const {return thetaCM;}    ///< scattering angle in the relative (CM) frame (rad)
  Double_t GetImpactParameter() const {return bimp;}   ///< effective \f$b=|\vec L|/(\mu v_\infty)\f$ (fm)
  /// @}

  /// \name Asymptotes (lab)
  /// Lab velocities of each particle at \f$t\to\mp\infty\f$ (c): the centre-of-mass
  /// velocity plus its share of the relative velocity. A fixed particle has
  /// the zero vector.
  /// @{
  TVector3 GetV1In()        const {return v1in;} ///< particle 1, incoming
  TVector3 GetV1Out()       const {return v1out;} ///< particle 1, outgoing
  TVector3 GetV2In()        const {return v2in;} ///< particle 2, incoming
  TVector3 GetV2Out()       const {return v2out;} ///< particle 2, outgoing
  /// @}
  /// \brief Lab scattering angle of particle 1 (rad).
  ///
  /// The angle between its own asymptotic incoming and outgoing lab velocity.
  /// If that particle is at rest before the collision (target at rest, fixed
  /// centre), it is measured from the incoming direction of the OTHER
  /// particle (recoil angle); 0 if both are 0. 0 for a fixed particle and
  /// when there is no asymptote.
  ///
  /// \warning Meaningless when that particle's outgoing speed is ~0 (equal
  /// masses head-on). "At rest" means \f$|v_{in}|\le 10^{-9}|v_{in}^{other}|\f$:
  /// for a target that is at rest at a FINITE distance, its asymptotic
  /// incoming velocity is not zero.
  Double_t GetThetaLab1()   const {return thetaLab1;}
  /// \brief Lab scattering angle of particle 2 (rad); same rules as GetThetaLab1().
  Double_t GetThetaLab2()   const {return thetaLab2;}

  /// \brief Text summary (status, kappa, mu, E, |L|, e, rmin, b, vinf, theta, outgoing velocities).
  /// \param mute if true, only build the string; otherwise also print it on \c std::cout
  /// \return the text
  TString  Print(Bool_t mute=kFALSE) const;

  /// \name Closed-form relations (for validation; angles in rad)
  /// Static helpers with the textbook Rutherford formulas, independent of any object.
  /// @{

  /// \brief Rutherford scattering angle in the CM frame.
  ///
  /// \f[ \theta = 2\arctan\frac{\kappa}{\mu v_\infty^2 b} \f]
  /// Signed: negative for an attractive \c kappa (deflection toward the
  /// centre); \f$|\theta|\f$ is the angle between incoming and outgoing
  /// directions. \f$b=0\f$ gives \f$\pm\pi\f$. The sign of \c b is ignored
  /// (\f$|b|\f$ is used).
  /// \param b impact parameter (fm)
  /// \param vinf relative speed at infinity (c)
  /// \param mu reduced mass (amu)
  /// \param kappa interaction constant
  static Double_t RutherfordAngle(Double_t b, Double_t vinf, Double_t mu, Double_t kappa);
  /// \brief Inverse of RutherfordAngle().
  ///
  /// \f[ b = \frac{|\kappa|}{\mu v_\infty^2\tan(|\theta|/2)} \f]
  /// Gives infinity for \f$\theta=0\f$.
  /// \param theta scattering angle (rad); the sign is ignored
  /// \param vinf relative speed at infinity (c)
  /// \param mu reduced mass (amu)
  /// \param kappa interaction constant
  static Double_t RutherfordImpactParameter(Double_t theta, Double_t vinf, Double_t mu, Double_t kappa);
  /// \brief Closest approach for impact parameter \c b.
  ///
  /// \f[ r_{min} = a + \sqrt{a^2 + b^2},\qquad a = \frac{\kappa}{\mu v_\infty^2} \f]
  /// (evaluated as \f$b^2/(\sqrt{a^2+b^2}-a)\f$ for \f$a<0\f$ to avoid cancellation).
  /// \param b impact parameter (fm)
  /// \param vinf relative speed at infinity (c)
  /// \param mu reduced mass (amu)
  /// \param kappa interaction constant
  static Double_t RutherfordClosestApproach(Double_t b, Double_t vinf, Double_t mu, Double_t kappa);
  /// \brief Head-on (\f$b=0\f$) closest approach, repulsive.
  ///
  /// \f[ d = \frac{\kappa}{E} = \frac{2\kappa}{\mu v_\infty^2} \f]
  /// Negative (meaningless) for an attractive \c kappa.
  /// \param vinf relative speed at infinity (c)
  /// \param mu reduced mass (amu)
  /// \param kappa interaction constant
  static Double_t HeadOnClosestApproach(Double_t vinf, Double_t mu, Double_t kappa);
  /// \brief Lab angle of the projectile (mass \c m1) on a target (mass \c m2) at rest, from the CM angle.
  ///
  /// \f[ \tan\theta_{lab} = \frac{\sin\theta}{\cos\theta + m_1/m_2} \f]
  /// \param thetaCM scattering angle in the CM frame (rad)
  /// \param m1 projectile mass (amu)
  /// \param m2 target mass (amu)
  static Double_t LabAngle(Double_t thetaCM, Double_t m1, Double_t m2);
  /// \brief Lab angle of the recoiling target at rest, \f$(\pi-\theta)/2\f$.
  /// \param thetaCM scattering angle in the CM frame (rad)
  static Double_t RecoilAngle(Double_t thetaCM);
  /// \brief Angle between two vectors, \f$\mathrm{atan2}(|\vec a\times\vec b|, \vec a\cdot\vec b)\f$, accurate at 0 and \f$\pi\f$.
  /// \return the angle in rad, 0 if one of the vectors is zero
  static Double_t Angle(const TVector3 & a, const TVector3 & b);
  /// @}

private:
  void Reset();     ///< set every result and input to its empty value (status kInvalid)
  void Solve();     ///< from r, v, mu, kappa (set by Compute)
  void SetLab();    ///< lab velocities and angles from the relative asymptotes

  status   st          = kInvalid; ///< orbit classification
  Bool_t   fixedCentre = kFALSE;   ///< exactly one body fixed
  Bool_t   fixed1      = kFALSE;   ///< particle 1 is a fixed centre
  Bool_t   fixed2      = kFALSE;   ///< particle 2 is a fixed centre
  Double_t m1          = 0.;       ///< mass of particle 1 (amu)
  Double_t m2          = 0.;       ///< mass of particle 2 (amu)
  Double_t mu          = 0.;       ///< reduced (or moving) mass (amu)
  Double_t kappa       = 0.;       ///< interaction constant
  TVector3 r, v, vcm;              ///< relative position, relative velocity, centre-of-mass velocity

  Double_t energy      = 0.;       ///< orbital energy
  TVector3 lvec, avec, evec, phat; ///< angular momentum, Laplace-Runge-Lenz, eccentricity vector, periapsis direction
  Double_t ecc         = 0.;       ///< eccentricity
  Double_t rmin        = 0.;       ///< closest approach

  Double_t vinf        = 0.;       ///< relative speed at infinity
  TVector3 vinRel, voutRel;        ///< incoming / outgoing relative velocity
  Double_t thetaCM     = 0.;       ///< CM scattering angle
  Double_t bimp        = 0.;       ///< effective impact parameter

  TVector3 v1in, v1out, v2in, v2out; ///< lab asymptotic velocities
  Double_t thetaLab1   = 0.;       ///< lab angle of particle 1
  Double_t thetaLab2   = 0.;       ///< lab angle of particle 2
};

#endif
