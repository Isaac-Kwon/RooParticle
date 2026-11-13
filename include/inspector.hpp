#ifndef H_INSPECTOR
#define H_INSPECTOR 1;

#include "TROOT.h"
#include "particle.hpp"
#include "event.hpp"

class event;

/// \brief Inspector : one terminating condition of event::DeriveInspect().
///
/// The event stops when its inspectors are satisfied. Normal inspectors are combined by
/// event::SetInspectMode(): AND by default (all must be satisfied), or OR (any). The AND
/// default is the design: for example "the particle reached the asymptote (angle between
/// dX and dV below a value) AND went through a minimum number of steps".
/// An inspector set as abort (SetAbort()) stops the event by itself, whatever the others say
/// (close encounter, numerical failure, ...).
///
/// Two kinds exist: inspectorP looks at two chosen particles, and inspectorE at the whole
/// event (all particles). In an event with many bodies one can add several inspectorP, for
/// different pairs, and inspectorE together; they are all combined by the inspect mode.
/// The inspector is owned by the event it was added to.
///
/// The method is a name ("DEG", ...) or the equal enumerator (methodID). The comparison
/// direction is in the name: a trailing G means greater than the value, L less than.
/// An unknown name is accepted when set but gives kInvalid, and event::DeriveInspect()
/// then refuses to start.
class inspector{
  public:
  /// \brief Method identifier.
  ///
  /// The string names (in the tables of inspectorP and inspectorE) and these enumerators are
  /// the same thing : "DEG" is kDEG, ... An inspectorP accepts the first group, an
  /// inspectorE the second. kInvalid is what an unknown name is parsed to.
  enum methodID {
    kInvalid = 0, ///< unknown method name
    // inspectorP
    kCOS,  ///< "COS": cosine of the angle between dX and dV < value
    kDEG,  ///< "DEG": angle between dX and dV in degree < value
    kRAD,  ///< "RAD": angle between dX and dV in radian < value
    kCNT,  ///< "CNT": number of derived steps > value
    kSPDG, ///< "SPDG": relative speed > value
    kSPDL, ///< "SPDL": relative speed < value
    kDSTG, ///< "DSTG": distance > value
    kDSTL, ///< "DSTL": distance < value
    kVRG,  ///< "VRG": radial velocity > value
    kVRL,  ///< "VRL": radial velocity < value
    // inspectorE
    kKEG,  ///< "KEG": net kinetic energy > value
    kKEL,  ///< "KEL": net kinetic energy < value
    kPEG,  ///< "PEG": net potential energy > value
    kPEL,  ///< "PEL": net potential energy < value
    kPKG,  ///< "PKG": abs(net PE / net KE) > value
    kPKL,  ///< "PKL": abs(net PE / net KE) < value
    kEDG,  ///< "EDG": relative drift of the net energy > value
    kEDL,  ///< "EDL": relative drift of the net energy < value
    kTMG,  ///< "TMG": evolved time > value
    kTML,  ///< "TML": evolved time < value
    kDXG,  ///< "DXG": largest step count > value
    kDXL,  ///< "DXL": smallest step count < value (Evaluate() gives the largest, see inspectorE)
    kDNG,  ///< "DNG": largest step count > value (Evaluate() gives the smallest, see inspectorE)
    kDNL,  ///< "DNL": smallest step count < value
    kDSG,  ///< "DSG": smallest distance between two particles > value
    kDSL,  ///< "DSL": smallest distance between two particles < value
    kACG,  ///< "ACG": largest acceleration > value
    kACL,  ///< "ACL": largest acceleration < value
    kVOUT, ///< "VOUT": any moving particle is outside the volume
    kVOUTA,///< "VOUTA": all moving particles are outside the volume
    kASYD, ///< "ASYD": remaining deflection of every moving particle < value (rad), all pairs receding
    kASYV  ///< "ASYV": remaining relative velocity change of every moving particle < value, all pairs receding
  };

  /// \brief Parse a method name.
  /// \param name e.g. "DEG" (exact, upper case)
  /// \return the identifier, kInvalid if unknown.
  static methodID ParseMethod(TString name);

  /// \brief Name of a method identifier.
  /// \return the name, "" for kInvalid.
  static TString  MethodName(methodID id);

  /// \brief Make an inspector without method (invalid until SetMethod()).
  inspector(){}
  virtual ~inspector(){}

  /// \brief Whether the condition is satisfied now (the base class : never).
  virtual Bool_t  Inspect(){return kFALSE;}

  /// \brief The quantity that is compared with the criterion (the base class : 0).
  virtual Double_t Evaluate(){return 0.;}

  /// \brief Whether the method name is defined for this kind of inspector (the base class : never).
  virtual Bool_t  IsValid(){return kFALSE;}

  /// \brief Set the criterion value to compare with (in the unit of the quantity).
  void    SetValue(Double_t val_){val=val_;}

  /// \brief Set the method by name.
  /// \param methodname e.g. "DEG"
  ///
  /// The name is parsed once here; Inspect()/Evaluate() dispatch on the enumerator.
  /// An unknown name is kept as it is (GetMethod() returns it) but gives kInvalid.
  void    SetMethod(TString methodname){ftn=methodname; mtd=ParseMethod(methodname);}

  /// \brief Set the method by C string (see SetMethod(TString)).
  void    SetMethod(const char * methodname){SetMethod(TString(methodname));}

  /// \brief Set the method by identifier.
  void    SetMethod(methodID id_){mtd=id_; ftn=MethodName(id_);}

  /// \brief The method name as given.
  TString GetMethod(){return ftn;}

  /// \brief The method identifier (kInvalid for an unknown name).
  methodID GetMethodID(){return mtd;}

  /// \brief The criterion value (default -1).
  Double_t GetCriterion(){return val;}

  /// \brief Make this inspector an abort inspector.
  /// \param abort_ kTRUE (default): when satisfied it stops the event by itself, whatever the others say.
  void    SetAbort(Bool_t abort_=kTRUE){abort=abort_;}

  /// \brief Whether this is an abort inspector.
  Bool_t  IsAbort(){return abort;}

  /// \brief Describe the inspector as text (the base class : empty).
  /// \param element kTRUE : include the particles (inspectorP)
  /// \param mechanics not used
  /// \param mute kFALSE (default) also prints to stdout
  /// \param pprint not used
  virtual TString Print(Bool_t element=kFALSE, Bool_t mechanics=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);

  /// \brief Whether this is an inspectorP.
  virtual Bool_t  IsInspectorP(){return kFALSE;}

  /// \brief Whether this is an inspectorE.
  virtual Bool_t  IsInspectorE(){return kFALSE;}

  protected:
  TString ftn="";             ///< method name as given
  methodID mtd=kInvalid;      ///< parsed method
  Double_t val=-1;            ///< criterion value
  Bool_t abort=kFALSE;        ///< abort inspector

  private:
  // inspectorT(particle *p1_, particle *p2_, TString method_, Double_t val_=-1){};
  // inspectorT(event* event_, TString method_, Double_t val_=-1){};

};

/// \brief Inspector on a pair of particles.
///
/// dX and dV are the relative position (fm) and velocity (c) of the two particles,
/// \f$\vec{dX} = \vec x_1 - \vec x_2\f$. With many bodies, each inspectorP looks only at its
/// own two particles; the other bodies are ignored.
///
/// | Method | Satisfied when |
/// |---|---|
/// | `"COS"` | cosine of the angle between dX and dV is less than the value |
/// | `"DEG"` | angle between dX and dV, in degree, is less than the value |
/// | `"RAD"` | angle between dX and dV, in radian, is less than the value |
/// | `"CNT"` | number of derived steps is greater than the value |
/// | `"SPDG"` | relative speed (norm of dV) is greater than the value |
/// | `"SPDL"` | relative speed (norm of dV) is less than the value |
/// | `"DSTG"` | distance (norm of dX) is greater than the value |
/// | `"DSTL"` | distance (norm of dX) is less than the value |
/// | `"VRG"` | radial velocity dX.dV/norm(dX) is greater than the value (value 0 : receding) |
/// | `"VRL"` | radial velocity dX.dV/norm(dX) is less than the value (value 0 : approaching) |
///
/// The enumerators are kCOS kDEG kRAD kCNT kSPDG kSPDL kDSTG kDSTL kVRG kVRL
/// (inspector::methodID), so `new inspectorP(p1, p2, inspector::kDEG, 1)` is the same as
/// `new inspectorP(p1, p2, "DEG", 1)`.
///
/// The angle is \f$\mathrm{atan2}(|\vec{dX}\times\vec{dV}|, \vec{dX}\cdot\vec{dV})\f$ : accurate
/// also near 0 and 180 degree (acos of the cosine is not). "COS" is clamped to [-1, 1].
/// A small "DEG" means dX and dV are parallel, i.e. the particles fly apart: it is the
/// asymptote condition after a scattering. "CNT" is the smaller of the step counts of the
/// two particles (an invincible particle counts steps too).
///
/// \par Zero vectors
/// If dX or dV is exactly the zero vector the angle / cosine (and the radial velocity for
/// dX = 0) is undefined: Evaluate() returns NaN and Inspect() is kFALSE, i.e. the inspector
/// is simply not satisfied (also for the "<" methods).
class inspectorP: public inspector{ //Calculator, Static
  public:
  /// \brief Make an inspector on two particles, method by name.
  /// \param p1_ first particle (of the event)
  /// \param p2_ second particle (of the event)
  /// \param method_ method name, see the table
  /// \param val_ criterion value (angle in degree for "DEG", radian for "RAD", distance in
  ///        fm, speed in c, count for "CNT")
  ///
  /// The particles must belong to the event the inspector is added to (event::event(event&)
  /// relinks it by index).
  inspectorP(particle *p1_, particle *p2_, TString method_, Double_t val_=-1);

  /// \brief Make an inspector on two particles, method by enumerator (see above).
  inspectorP(particle *p1_, particle *p2_, methodID method_, Double_t val_=-1);

  /// \brief Whether the condition is satisfied now.
  /// \return the comparison of Evaluate() with the criterion, in the direction given by the
  ///         table; kFALSE for an invalid method or a zero vector.
  Bool_t  Inspect();

  /// \brief The compared quantity: cosine, angle (degree or radian), step count, relative
  ///        speed (c), distance (fm) or radial velocity (c). NaN for a zero vector (angle,
  ///        cosine, radial velocity); for an invalid method it prints a message and gives 0.
  Double_t Evaluate();

  /// \brief Whether the method is one of the inspectorP methods.
  Bool_t  IsValid();

  /// \brief Set both particles.
  void SetParticles(particle * p1_, particle * p2_){p1 = p1_; p2 = p2_;}
  /// \brief Set the first particle.
  void SetParticle1(particle * p1_){p1 = p1_;}
  /// \brief Set the second particle.
  void SetParticle2(particle * p2_){p2 = p2_;}

  /// \brief The first particle.
  particle* GetParticle1(){return p1;}
  /// \brief The second particle.
  particle* GetParticle2(){return p2;}

  /// \brief Describe as text "[ Inspector Particle: MTD, CRT (criterion), VAL (value), RST (result) ]".
  /// \param particles kTRUE also prints the two particles
  /// \param mechanics not used
  /// \param mute kFALSE (default) also prints to stdout
  /// \param pprint not used
  TString Print(Bool_t particles=kFALSE, Bool_t mechanics=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);

  // Bool_t SetInitial(Bool_t index=0, TString mtd="", Double_t var=0.1);

  /// \brief Move one particle back along its velocity so that the present condition is
  ///        the one of the given angle (setting the initial position from the present state).
  /// \param index false : move p2 ; true : move p1
  /// \param mtd "COS", "DEG" or "RAD" (empty : the method of the inspector)
  /// \param var the target angle margin
  /// \return kTRUE, or kFALSE for another method.
  /// \warning Unfinished and not trustworthy: the function itself prints "needs to be fixed.
  ///          do not use without any instruction" on stderr. It changes the position of a
  ///          particle.
  /// \todo Fix or remove (see TODO.md).
  Bool_t SetInitial_D(Bool_t index=0, TString mtd="", Double_t var=0.1);

  /// \brief Always kTRUE.
  Bool_t IsInspectorP(){return kTRUE;}

  protected:
  /// \brief Whether dX (or, if needV, dV) is the zero vector.
  Bool_t   IsDegenerate(Bool_t needV=kTRUE);
  /// \brief Cosine of the angle between dX and dV, clamped to [-1, 1] (NaN for a zero vector).
  Double_t CalculateCosine();
  /// \brief Angle between dX and dV in degree (radian if rad), NaN for a zero vector.
  Double_t CalculateAngle(Bool_t rad=kFALSE);
  /// \brief Smaller of the step counts of the two particles.
  Long_t  CountDeriving();
  /// \brief Relative speed, norm of dV (c).
  Double_t CaculateSPD();
  /// \brief Distance, norm of dX (fm).
  Double_t CalculateDistance();
  /// \brief dX.dV/norm(dX) (c), positive when receding; NaN if dX = 0.
  Double_t CalculateRadialVelocity();

  private:
  particle *p1, *p2; ///< the two particles (not owned)
};


/// \brief Inspector on the whole event.
///
/// All particles of the event are taken into account, so these methods suit an event with
/// many bodies. "Moving" below means non-invincible.
///
/// | Method | Satisfied when |
/// |---|---|
/// | `"KEG"` / `"KEL"` | net kinetic energy of the moving particles (event::GetNetKE(), non-relativistic) is greater / less than the value |
/// | `"PEG"` / `"PEL"` | net potential energy (event::GetNetPE()) is greater / less than the value |
/// | `"PKG"` / `"PKL"` | abs(net PE / net KE) is greater / less than the value (PKL : the interaction became negligible); with KE = 0 the ratio is the largest double |
/// | `"EDG"` / `"EDL"` | relative drift of the net energy abs(E-E0)/abs(E0) is greater / less than the value (E0 : energy at the first evaluation; EDG as abort : numerical failure) |
/// | `"TMG"` / `"TML"` | evolved time (event::GetTime(), fm/c) is greater / less than the value |
/// | `"DXG"` | largest step count (event::GetDeriveN(kFALSE)) is greater than the value |
/// | `"DXL"` | smallest step count (event::GetDeriveN(kTRUE)) is less than the value |
/// | `"DNG"` | largest step count is greater than the value |
/// | `"DNL"` | smallest step count is less than the value |
/// | `"DSG"` / `"DSL"` | smallest distance between any two particles (fm) is greater / less than the value (DSL : close encounter) |
/// | `"ACG"` / `"ACL"` | largest acceleration of the moving particles is greater / less than the value (ACL : forces became negligible) |
/// | `"VOUT"` | any moving particle is outside the event volume (value not used; never satisfied without a volume) |
/// | `"VOUTA"` | all moving particles are outside the event volume (value not used) |
/// | `"ASYD"` | the deflection every moving particle still has ahead of it, from now to infinity, is less than the value (rad), and every interacting pair is receding |
/// | `"ASYV"` | the same for the relative velocity change max\|v_inf - v\|/\|v_inf\| |
///
/// The enumerators are kKEG kKEL kPEG kPEL kPKG kPKL kEDG kEDL kTMG kTML kDXG kDXL kDNG
/// kDNL kDSG kDSL kACG kACL kVOUT kVOUTA kASYD kASYV (inspector::methodID).
///
/// \par Asymptote reached, for any number of bodies
/// "ASYD" and "ASYV" are the terminating conditions meant for three and more bodies: they
/// ask directly how far the outgoing velocities still are from their values at infinity,
/// through the first-order estimate asymptoteN::RemainingDeflection() /
/// asymptoteN::RemainingChange(). They are never satisfied while any interacting pair is
/// still approaching. One evaluation costs about as much as a few force sweeps, so inspect
/// every 10 to 100 steps (iperiod of event::DeriveInspect()).
///
/// \par Evaluate() against Inspect()
/// Evaluate() returns the number compared in the table, with one asymmetry: for the step
/// count, Inspect() uses the largest count for "DXG" and "DNG" and the smallest for "DXL" and
/// "DNL", while Evaluate() returns the largest count for "DXG"/"DXL" and the smallest for
/// "DNG"/"DNL". So the value printed for "DXL" or "DNG" is not the one Inspect() compared.
/// For "VOUT"/"VOUTA" Evaluate() returns 1 or 0 (the result itself).
///
/// \par Energy drift
/// "EDG"/"EDL" use the relativistic energy event::GetNetERel() when the relativistic mode is
/// really in effect (event::SetRelativistic() with an integrator other than kTaylor), and
/// event::GetNetE() otherwise. The reference E0 is taken at the first evaluation after
/// the inspector was made, set to another event, or reset with ResetReference(); with
/// E0 = 0 the absolute value of the energy is used. "KEG"/"KEL"/"PKG"/"PKL" always use the
/// non-relativistic kinetic energy.
///
/// \par Copies
/// When an event is copied, its inspectorE are remade to watch the new event, and E0 is
/// forgotten.
class inspectorE: public inspector{
  public:
  /// \brief Make an inspector on an event, method by name.
  /// \param event_ the event to watch (normally the event the inspector is added to)
  /// \param method_ method name, see the table
  /// \param val_ criterion value (energies in amu c^2, time in fm/c, distance in fm, acceleration in c^2/fm, ratio and drift without unit)
  inspectorE(event* event_, TString method_, Double_t val_=-1);

  /// \brief Make an inspector on an event, method by enumerator.
  inspectorE(event* event_, methodID method_, Double_t val_=-1);

  /// \brief Set the event to watch (forgets the reference energy E0).
  void     SetEvent(event* event_){evt=event_; hasRef=kFALSE;}

  /// \brief The watched event.
  event *  GetEvent(){return evt;}

  /// \brief Whether the condition is satisfied now (kFALSE for an invalid method, with a message on stderr).
  Bool_t   Inspect();

  /// \brief The compared quantity, see the table (the first "EDG"/"EDL" call sets E0 and gives 0).
  Double_t Evaluate();

  /// \brief Whether the method is one of the inspectorE methods.
  Bool_t   IsValid();

  /// \brief Forget E0 of "EDG"/"EDL"; the next evaluation takes the present energy as E0.
  void     ResetReference(){hasRef=kFALSE;}

  /// \brief Describe as text "[ Inspector Event: MTD, CRT, VAL, RST (, ABORT) ]".
  /// \param element not used
  /// \param mechanics not used
  /// \param mute kFALSE (default) also prints to stdout
  /// \param pprint not used
  TString  Print(Bool_t element=kFALSE, Bool_t mechanics=kFALSE, Bool_t mute=kFALSE, Bool_t pprint=kFALSE);

  /// \brief Always kTRUE.
  Bool_t   IsInspectorE(){return kTRUE;}

  protected:

  private:
  event *evt;              ///< the watched event (not owned)
  Double_t ref = 0.;       ///< reference energy E0 of "EDG"/"EDL"
  Bool_t   hasRef = kFALSE; ///< E0 is set
  // TString ftn=TString("");
  // Double_t val=-1;

};

#endif
