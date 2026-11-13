#ifndef H_ASYMPTOTEN
#define H_ASYMPTOTEN 1

#include "vector"

#include "TROOT.h"
#include "TString.h"
#include "TVector3.h"

class event;

////////////////////////////////////////////////////////////////////////////
//
// asymptoteN : asymptotic (t -> +-infinity) velocities of an N-body event,
//              with an error estimate
//
// A run stops at a finite distance, where the interaction is not negligible.
// For two bodies with a pure 1/r potential the class asymptote is exact. For
// N bodies there is no closed form, so asymptoteN returns, for every moving
// particle, the asymptotic velocity v_inf AND an estimate of its error, for
// any number of bodies and any of the pair forces of force.hpp. The caller's
// event is never modified (the integrating methods work on a private copy).
//
// What "asymptotic state" means here : VELOCITIES ONLY. With 1/r forces the
// positions never settle on straight lines (r = v t + (kappa/mu v^2) ln t + ...,
// a logarithmic phase), while the velocities converge like 1/t. Two bodies
// that leave bound to each other have no individual asymptotic velocity : they
// are reported as a cluster with its centre-of-mass velocity (see Clusters).
//
// Methods (Compute(ev, method, direction))
//   kBorn      First order (straight-line impulse). Every particle moves on a
//              straight line x + v t; the impulse of each pair force is
//              integrated to t = infinity : closed form for unsoftened Coulomb
//              and gravity, I = kappa/(w R) [ b_vec/(R + s) + w_hat ]
//              (R = |r|, s = r.w_hat, b_vec = r - s w_hat, w = |v1 - v2|),
//              adaptive Gauss-Legendre quadrature of force::Evaluate() along
//              the line for the other forces. O(N^2) and cheap : this is the
//              residual / termination criterion. Works in the relativistic
//              mode (it integrates dp/dt along the lab-frame straight lines).
//              Error : O(eta^2) per pair, eta ~ |U|/(mu w^2) ~ 1/R ; the
//              measured order of the velocity error in 1/R is 2.0.
//   kPairwise  (default) Every pair is solved as an isolated two-body problem
//              EXACTLY, and the velocity changes of all pairs are added :
//              pure 1/r pairs with the class asymptote, any other central
//              potential (Coulomb + Woods-Saxon, Yukawa, softened forces) with
//              the orbit integral  dphi = int h / (r^2 sqrt(2(E-U)/mu - h^2/r^2)) dr
//              (adaptive quadrature, turning point included, so an approaching
//              pair, even one passing through a nucleus, is fine too).
//              Exact for two bodies (1e-16 closed form, 1e-15 quadrature);
//              for N >= 3 the error is the three-body cross term (pair ik
//              perturbed by the future kicks of the others), also O(1/R^2),
//              with a 10 - 60 times smaller coefficient than kBorn in the tests.
//              Momentum is conserved exactly ; the energy defect at infinity
//              is O(cross term) (GetEnergyInf() - GetEnergyNow()).
//              A receding pair whose non-1/r forces are below 1e-17 of its
//              potential and force (outside a nuclear range) is solved in closed form.
//   kAuto      kPairwise on a private copy of the event that is continued with
//              the adaptive RK45 in doubling time blocks until the estimated
//              error is below SetAutoTolerance(). 1 - 11 blocks in the tests ;
//              at the default 1e-9 it costs about as much as kBrute, it pays
//              off for looser tolerances.
//   kBrute     Reference : a private copy is continued with RK45 in doubling
//              time blocks t_k ; the velocities are extrapolated with the 1/t
//              law (Richardson, v_inf = (t2 v2 - t1 v1)/(t2 - t1)) until two
//              successive extrapolations agree within SetBruteTolerance().
//              Its residual error is O(ln t / t^2). Independent of the
//              analytic methods (it never uses them), slow compared to them.
//
// Direction
//   kOutgoing  t -> +infinity (default).
//   kIncoming  t -> -infinity : computed by time reversal (all velocities of
//              a copy flipped, outgoing asymptote, flipped back). Valid for
//              velocity-independent forces, which all forces of force.hpp are.
//              It gives the TRUE incoming velocities (and so the effective
//              impact parameters) of a run started at a finite distance.
//
// Error estimate (GetError(i), in c ; a 1-sigma-like magnitude, not a strict bound)
//   kPairwise  err_a = C sum_b (mu_ab/m_a) eta_ab |dw_ab^ext|            (cross term)
//                    + quadrature error (two orders compared)
//                    + multipole error of clusters (members vs point composite)
//              eta_ab = max(|U|, r |F_r|) / (2 E_ab) at r = current distance
//              (at the closest approach still ahead for an approaching pair),
//              dw_ab^ext = change of the relative velocity of a,b caused by the
//              OTHER bodies (from the same pairwise solution), C = SetCrossFactor()
//              (default 2 : one eta from the displaced relative orbit, one from
//              the shifted relative velocity it is solved with). Zero for two
//              bodies, where only the quadrature and round-off (4e-15 of the
//              largest velocity) terms remain.
//   kBorn      the kPairwise estimate + |dv_Born - dv_Pairwise| (second-order
//              size of the pair terms : for two bodies it IS the error).
//   kAuto      the kPairwise estimate at the last block + integration error.
//   kBrute     |difference of the last two extrapolations| + integration error.
//   Integration error of a private copy : (rel |v| + abs) sqrt(n substeps) with
//   the RK45 tolerances of SetIntegrationTolerance().
//   Measured by main_asymptoteN_check.cpp, actual error / estimate, 3 and 4
//   bodies (fixed and free targets, Coulomb and Coulomb + Woods-Saxon, stops
//   from 15 fm to 23000 fm, outgoing and incoming) : kPairwise and kAuto
//   <= 0.23 (the estimate is 4 - 9 times too large), kBorn <= 0.99 (tight :
//   it is dominated by the exactly known pair term), kBrute <= 0.34.
//
// Residual and termination criterion
//   GetResidual(i) = v_inf - v_now, GetDeflection(i) = remaining deflection
//   angle (rad) of particle i, GetMaxDeflection() over the moving particles.
//   RemainingDeflection(ev) is the cheap static form (kBorn, no estimate) for
//   a stop condition : "stop when every particle will still turn by < X rad"
//   (use together with IsSeparating() : all interacting pairs receding).
//
// Clusters (bound subsystems)
//   Bodies are merged while some pair of bodies has a negative pair energy
//   E_ab = mu_ab |V_a - V_b|^2 / 2 + sum U_ij (members) in its own CM frame.
//   A cluster is then treated as ONE body (total mass and charge at its centre
//   of mass) and every member gets the cluster velocity, with flag kInCluster.
//   Its multipole interaction (e.g. the dipole of a neutral bound pair) is not
//   in v_inf : it enters the error estimate as |first-order impulse of the
//   members, moving rigidly - that of the point composite|.
//   GetClusterEnergy(c) = internal energy (kinetic in the CM frame + internal
//   potential), GetClusterIsolation(c) = |tidal / multipole coupling to the
//   other bodies| / |internal energy|; above SetIsolation() (default 1e-2) :
//   kNotIsolated, the binding may still be broken (or made) later.
//   Bound to an invincible particle : kBoundToFixed, v_inf = 0 for its members.
//
// Flags (per particle GetFlags(i), OR of all GetFlags())
//   kApproaching   an interacting pair of it is still approaching : handled,
//                  but the cross-term estimate is less certain.
//   kStrong        some pair has C eta > 0.5 and is perturbed by other bodies :
//                  the perturbative estimate is NOT trustworthy (stop later,
//                  or use kAuto / kBrute).
//   kInCluster     member of a bound cluster : v_inf is the cluster CM velocity.
//   kNotIsolated   that cluster is not isolated (see above).
//   kBoundToFixed  bound to an invincible particle : no asymptote, v_inf = 0.
//   kRelativistic  relativistic event : kBorn (first order) is used instead of
//                  kPairwise, which is non-relativistic.
//   kNoPairwise    a force is not central or the event is not pair symmetric :
//                  first order (kBorn) used.
//   kMarginal      a pair with E ~ 0 (v_inf ~ 0) : slow convergence, unreliable.
//   kSingular      radial fall into an attractive singular centre, or a
//                  first-order impulse that diverges (head-on approach).
//   kNotConverged  kAuto / kBrute reached SetMaxBlocks() before the tolerance.
//   kInvalid       no event, no moving particle, zero mass ...
//   IsReliable() : none of kStrong, kNotIsolated, kMarginal, kSingular,
//   kNotConverged, kInvalid.
//
// Validity and failure modes
//   * Forces must be pairwise, velocity independent and vanish at infinity.
//     kPairwise needs them central (all forces of force.hpp are); a
//     non-central force falls back to kBorn (flag kNoPairwise).
//   * Particles are evaluated through stand-ins carrying mass, charge and
//     position (EMparticle) : a force reading other particle data is not
//     supported (the event copy constructor makes the same assumption).
//   * Not trustworthy while particles still approach each other strongly,
//     for slow relative motion (eta not small) and for clusters that are not
//     isolated : the flags say so. kAuto handles the first two by integrating on.
//   * Relativistic mode : first order (kBorn) only, error O(eta^2) estimated
//     from the non-relativistic second-order size; kBrute integrates the
//     relativistic copy and is exact up to its tolerance. Retardation is
//     ignored (as in the event).
//   * oddevent (sequential update) is not modelled : the copies are event.
//   * Invincible particles are fixed (their stored velocity is ignored, as the
//     event never moves them) and are excluded from momentum and energy sums.
//
// Cost (2 - 4 bodies, measured by main_asymptoteN_check, one core)
//   kBorn / kPairwise, Coulomb            4 - 15 microseconds
//   kBorn / kPairwise, inside a nuclear range (quadrature)   0.2 - 0.4 ms
//   kBrute, unbound                       100 - 500 RK45 substeps, 0.3 - 1.3 ms
//                                         (to t ~ 1e9 fm/c ; tolerance 1e-13)
//   kBrute with a bound pair              its orbit is integrated : ~ 1e5 substeps
//   kAuto                                 0.3 - 3.5 ms at the default tolerance
//   The residual of the fixed-target three-body run : 8 microseconds per call.
//
// Example
//   asymptoteN asy(ev);                          // kPairwise, outgoing
//   for(Int_t i=0; i<ev->getNParticle(); i++){
//     TVector3 vout = asy.GetVInf(i);            // asymptotic velocity
//     Double_t dv   = asy.GetError(i);           // its error estimate (c)
//   }
//   if(!asy.IsReliable()) asy.Compute(ev, asymptoteN::kAuto);
//   asymptoteN in(ev, asymptoteN::kPairwise, asymptoteN::kIncoming);
//   // stop condition inside a loop :
//   Bool_t sep; if(asymptoteN::RemainingDeflection(ev, &sep) < 1e-6 && sep) break;
//
////////////////////////////////////////////////////////////////////////////
class asymptoteN{
public:
  enum method    { kBorn = 0, kPairwise, kAuto, kBrute };
  enum direction { kOutgoing = 0, kIncoming };
  enum flagBit {
    kApproaching  = 1<<0,
    kStrong       = 1<<1,
    kInCluster    = 1<<2,
    kNotIsolated  = 1<<3,
    kBoundToFixed = 1<<4,
    kRelativistic = 1<<5,
    kNoPairwise   = 1<<6,
    kMarginal     = 1<<7,
    kSingular     = 1<<8,
    kNotConverged = 1<<9,
    kInvalid      = 1<<10
  };
  static const UInt_t kUnreliable = kStrong | kNotIsolated | kMarginal | kSingular | kNotConverged | kInvalid;

  asymptoteN(){}
  asymptoteN(event* ev, method m = kPairwise, direction d = kOutgoing){Compute(ev, m, d);}

  // Compute from the present state of ev (ev is not changed). Returns kFALSE
  // when there is no result at all (kInvalid); check IsReliable() as well.
  Bool_t Compute(event* ev, method m = kPairwise, direction d = kOutgoing);

  // --- settings (keep their values across Compute calls)
  void SetAutoTolerance(Double_t tol){autoTol = tol;}         // kAuto : target of GetMaxRelError() (default 1e-9)
  void SetBruteTolerance(Double_t tol){bruteTol = tol;}       // kBrute : relative agreement of two extrapolations (1e-11)
  void SetIntegrationTolerance(Double_t rel, Double_t abs){rkRel = rel; rkAbs = abs;} // RK45 of the copies (1e-13, 1e-15)
  void SetMaxBlocks(Int_t n){maxBlocks = n;}                  // kAuto / kBrute : doubling blocks (default 60)
  void SetQuadratureOrder(Int_t n){nquad = (n<4) ? 4 : n;}    // Gauss points per panel (default 12)
  void SetIsolation(Double_t eps){isoThr = eps;}              // cluster isolation threshold (default 1e-2)
  void SetCrossFactor(Double_t c){crossC = c;}                // C of the cross-term estimate (default 2)
  void SetAnalytic(Bool_t on = kTRUE){analytic = on;}         // closed forms for pure 1/r pairs (off : quadrature everywhere)
  void SetErrorEstimate(Bool_t on = kTRUE){estimate = on;}    // kBorn : also solve the pairs for the estimate

  // --- results per particle (index of the event; invincible : v_inf = 0, error 0)
  Int_t    GetNParticle()          const {return (Int_t) vInf.size();}
  TVector3 GetV(Int_t i)           const {return vNow.at(i);}          // velocity of the given state (0 if fixed)
  TVector3 GetVBody(Int_t i)       const {return vBody.at(i);}         // GetV, or the cluster CM velocity for members
  TVector3 GetVInf(Int_t i)        const {return vInf.at(i);}          // asymptotic velocity (cluster CM velocity for members)
  Double_t GetError(Int_t i)       const {return err.at(i);}           // error estimate of GetVInf (c)
  TVector3 GetResidual(Int_t i)    const {return vInf.at(i) - vBody.at(i);} // velocity change still to come (or since t = -inf)
  Double_t GetDeflection(Int_t i)  const;                              // angle between GetVBody and GetVInf (rad)
  Double_t GetAngularError(Int_t i) const;                             // GetError / |GetVInf| (rad)
  UInt_t   GetFlags(Int_t i)       const {return flags.at(i);}
  Int_t    GetCluster(Int_t i)     const {return cluster.at(i);}       // -1 : not in a cluster
  Bool_t   IsFixed(Int_t i)        const {return fixed.at(i);}

  // --- results, whole event
  method    GetMethod()       const {return meth;}
  direction GetDirection()    const {return dir;}
  UInt_t   GetFlags()         const {return allFlags;}
  TString  GetFlagNames(UInt_t f) const;
  Bool_t   IsValid()          const {return !(allFlags & kInvalid);}
  Bool_t   IsReliable()       const {return !(allFlags & kUnreliable);}
  Bool_t   IsSeparating()     const {return separating;}   // every interacting pair receding (or bound in a cluster)
  Double_t GetMaxDeflection() const {return maxDefl;}      // largest remaining deflection (rad)
  Double_t GetMaxError()      const {return maxErr;}       // largest GetError (c)
  Double_t GetVelocityScale() const {return vScale;}       // largest |v_inf| (c)
  Double_t GetMaxRelError()   const {return (vScale>0.) ? maxErr/vScale : 0.;}
  Double_t GetMaxEta()        const {return maxEta;}       // largest pair strength eta = max(|U|, r|F|)/(2E)
  // Total energy now (moving particles' KE + all pair PE but fixed-fixed), and
  // the asymptotic one : sum of KE of the bodies at infinity + internal cluster
  // energies. Their difference measures the non-conservation of the method
  // (zero for two bodies). Non-relativistic kinetic energies.
  Double_t GetEnergyNow()     const {return eNow;}
  Double_t GetEnergyInf()     const {return eInf;}
  TVector3 GetMomentumInf()   const;                       // sum m v_inf over the moving particles

  // --- clusters
  Int_t    GetNClusters()     const {return (Int_t) clMembers.size();}
  std::vector<Int_t> GetClusterMembers(Int_t c) const {return clMembers.at(c);}
  Double_t GetClusterEnergy(Int_t c)    const {return clEnergy.at(c);}
  Double_t GetClusterIsolation(Int_t c) const {return clIso.at(c);}
  TVector3 GetClusterVInf(Int_t c)      const {return clVInf.at(c);}

  // --- cost of the last Compute
  Double_t GetRunTime()       const {return runTime;}      // wall-clock seconds
  Double_t GetIntegratedTime() const {return tInt;}        // fm/c continued by kAuto / kBrute
  Int_t    GetNBlocks()       const {return nBlocks;}
  Long_t   GetNSubsteps()     const {return nSub;}         // RK45 substeps (accepted + rejected)
  Long_t   GetNForceSweeps()  const {return 6*nSub + nBlocks;} // full force evaluations of the event

  TString  Print(Bool_t mute = kFALSE) const;

  // Cheap termination criterion : largest remaining deflection (rad) of a
  // moving particle (or cluster), first order, no error estimate.
  // separating (if given) : every interacting pair is receding.
  static Double_t RemainingDeflection(event* ev, Bool_t* separating = nullptr);
  // Same, largest |v_inf - v| / |v_inf| (relative velocity change still to come).
  static Double_t RemainingChange(event* ev, Bool_t* separating = nullptr);

  // the internal analysis (public for the implementation helpers only)
  struct pstate { Double_t m, q; TVector3 x, v; Bool_t fixed; };

private:
  void Reset(Int_t n);
  void Analyze(event* ev, const std::vector<pstate> & ps, method m, Bool_t clusterOnly = kFALSE);
  void RunIntegrating(event* ev, method m);
  void Finish(event* ev);

  // settings
  Double_t autoTol   = 1e-9;
  Double_t bruteTol  = 1e-11;
  Double_t rkRel     = 1e-13;
  Double_t rkAbs     = 1e-15;
  Int_t    maxBlocks = 60;
  Int_t    nquad     = 12;
  Double_t isoThr    = 1e-2;
  Double_t crossC    = 2.;
  Bool_t   analytic  = kTRUE;
  Bool_t   estimate  = kTRUE;

  // results
  method    meth = kPairwise;
  direction dir  = kOutgoing;
  std::vector<TVector3> vNow, vBody, vInf;
  std::vector<Double_t> err, mass;
  std::vector<UInt_t>   flags;
  std::vector<Int_t>    cluster;
  std::vector<Bool_t>   fixed;
  std::vector<std::vector<Int_t> > clMembers;
  std::vector<Double_t> clEnergy, clIso;
  std::vector<TVector3> clVInf;
  UInt_t   allFlags   = kInvalid;
  Bool_t   separating = kFALSE;
  Double_t maxDefl = 0., maxErr = 0., vScale = 0., maxEta = 0., etaApp = 0.;
  Double_t eNow = 0., eInf = 0.;
  Double_t runTime = 0., tInt = 0.;
  Int_t    nBlocks = 0;
  Long_t   nSub = 0;
};

#endif
