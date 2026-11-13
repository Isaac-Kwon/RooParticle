#ifndef H_FORCE
#define H_FORCE 1;

#include "TROOT.h"
#include "TVector3.h"
#include "particle.hpp"

/// \brief Base class of the pairwise forces.
///
/// `Evaluate(p1, p2)` is the force on p1 by p2 in amu c^2/fm and `Potential(p1, p2)` the
/// potential energy of the pair in amu c^2. Every subclass keeps
/// \f$\vec F_1 = -\nabla_{x_1} U\f$, so that the energy monitors (event::GetPE(),
/// the "EDG" inspector, ...) are consistent with the dynamics.
///
/// A force added to an event with event::AddForce() acts on every pair of particles of the
/// event, whatever their number: with N bodies the force on a particle is the sum of the
/// force by each of the other N-1 bodies, and of every force added. The force does not
/// know which pair it is asked about, so it must depend only on the two particles
/// (position, mass, charge). The event may call `Evaluate(p1, p2)` once per pair
/// (in the order of addition) and give the opposite to p2 (event::SetPairSymmetric()), so a force which is not
/// antisymmetric must turn that off. Forces are not owned by the event.
///
/// To make a new force, inherit it and define both functions.
class force{
public:
  /// \brief Virtual : a force deleted through force* must run the subclass destructor.
  virtual ~force(){}

  /// \brief Force on the first particle by the second.
  /// \param p1_ particle which feels the force
  /// \param p2_ particle which exerts the force
  /// \return force vector in amu c^2/fm. The base class prints "NULL Force" and returns the zero vector.
  virtual TVector3 Evaluate(particle* p1_, particle* p2_){
    TVector3 a_;
    std::cout<<"NULL Force"<<endl; return a_;
  }

  /// \brief Potential energy of the pair (amu c^2).
  /// \param p1_ first particle
  /// \param p2_ second particle
  /// \return the energy. The base class prints "NULL Force, Potential" and returns 0.
  virtual Double_t Potential(particle* p1_, particle* p2_){
    Double_t a_ = 0.;
    std::cout<<"NULL Force, Potential"<<endl; return a_;
  };
protected:
  /// \brief Only subclasses can be made.
  force(){}
private:
};

/// \brief Electric (Coulomb) force between the charges of two particles.
///
/// \f[ U = \frac{k q_1 q_2}{r} \f]
/// with the charges q1, q2 in e (particle::GetQ()) and
/// \f$k = 0.00154629\f$ amu c^2 fm / e^2 (about 1.44 MeV fm). Like charges repel, opposite
/// charges attract; a particle with charge 0 feels nothing. The force is
/// \f$\vec F_1 = k q_1 q_2 (\vec x_1 - \vec x_2)/r^3\f$, pointing from p2 to p1 when repulsive.
///
/// With softening eps [fm]:
/// \f[ U = \frac{k q_1 q_2}{\sqrt{r^2 + \epsilon^2}} \f]
/// and the force is its exact gradient. It keeps the force finite at r -> 0. eps = 0
/// (default) runs the original unsoftened code, bit for bit, and diverges (inf or NaN) when
/// two particles coincide. With many charged bodies the force on each is the vector sum of the
/// Coulomb forces by all the others.
class coulombForce: public force{
public:
  /// \brief Make the force with k = 0.00154629 and no softening.
  coulombForce():force(){}
  /// \brief Coulomb force on p1_ by p2_ (amu c^2/fm), see the class description.
  TVector3 Evaluate(particle* p1_, particle* p2_);
  /// \brief Coulomb potential energy of the pair (amu c^2), see the class description.
  Double_t Potential(particle* p1_, particle* p2_);

  /// \brief Set k in amu c^2 fm / e^2.
  void     SetK(Double_t k_){coulombK = k_;}
  /// \brief k in amu c^2 fm / e^2.
  Double_t GetK(){return coulombK;}
  /// \brief Set the softening length eps in fm (0 : none).
  void     SetSoftening(Double_t eps_){softening = eps_;}
  /// \brief The softening length eps in fm.
  Double_t GetSoftening(){return softening;}
protected:
  Double_t coulombK  = 0.00154629; ///< k in amu c^2 fm / e^2
  Double_t softening = 0.;         ///< softening length in fm
private:
};

/// \brief Newtonian gravity between the rest masses of two particles. Always attractive.
///
/// \f[ U = -\frac{G m_1 m_2}{r} \f]
/// with the rest masses m1, m2 in amu (also in the relativistic mode) and G in
/// amu^-1 c^2 fm. The default is the SI value converted to these units,
/// \f$G/c^2 \times 1\,\mathrm{amu} = 1.2331\times10^{-54}\,\mathrm{m} = 1.2331\times10^{-39}\,\mathrm{fm}\f$.
/// The force on p1 is \f$-G m_1 m_2 (\vec x_1 - \vec x_2)/r^3\f$, towards p2.
///
/// With softening eps [fm], r^2 is replaced by r^2 + eps^2 in the potential and the force is
/// its exact gradient, as in coulombForce. Unlike coulombForce, the softened formula is
/// always used, so eps = 0 gives inf or NaN for coinciding particles. With many bodies every
/// pair attracts, so the force on a body is the vector sum over all the others.
///
/// \note Until 2026-11 the default was 1.2331e-17, a wrong exponent (too strong by 1e22).
///       At nuclear scales gravity is about 1e-36 of the Coulomb force between two protons.
class universalForce: public force{ //Gravity
public:
  /// \brief Make the force with G = 1.2331e-39 and no softening.
  universalForce():force(){}
  /// \brief Gravitational force on p1_ by p2_ (amu c^2/fm).
  TVector3 Evaluate(particle* p1_, particle* p2_);
  /// \brief Gravitational potential energy of the pair (amu c^2).
  Double_t Potential(particle* p1_, particle* p2_);

  /// \brief Set the gravitational constant G.
  void     SetG(Double_t g_){gravityG = g_;}
  /// \brief The gravitational constant G.
  Double_t GetG(){return gravityG;}
  /// \brief Set the softening length eps in fm (0 : none).
  void     SetSoftening(Double_t eps_){softening = eps_;}
  /// \brief The softening length eps in fm.
  Double_t GetSoftening(){return softening;}
protected:
  Double_t gravityG  = 1.2331e-39; ///< gravitational constant (amu^-1 c^2 fm)
  Double_t softening = 0.;         ///< softening length in fm
private:
};

/// \brief Yukawa (meson exchange) potential, central, independent of charge and mass.
///
/// \f[ U = -\frac{g\, e^{-r/a}}{r} \f]
///  - g [amu c^2 fm] : coupling. g > 0 is attractive (nuclear), g < 0 repulsive.
///  - a [fm] : range. The default is the pion Compton wavelength
///    \f$\hbar/(m_\pi c) = 197.327/139.57 = 1.4138\f$ fm.
///
/// The default g = 0.1 amu c^2 fm (= 93.1 MeV fm) is only an illustrative strength, not a
/// fitted value: set it for the system. Like Coulomb, the force diverges at r -> 0 (inf or
/// NaN for coinciding particles). The force on p1 is along \f$\vec x_2 - \vec x_1\f$ for g > 0.
/// It acts on every pair of particles of the event (one particle = one body, not a sum over
/// nucleons), so with many bodies each one feels all the others.
class yukawaForce: public force{
public:
  /// \brief Make the force.
  /// \param g_ coupling in amu c^2 fm (default 0.1, illustrative)
  /// \param a_ range in fm (default 1.4138)
  yukawaForce(Double_t g_ = 0.1, Double_t a_ = 1.4138):force(){g = g_; a = a_;}
  /// \brief Yukawa force on p1_ by p2_ (amu c^2/fm).
  TVector3 Evaluate(particle* p1_, particle* p2_);
  /// \brief Yukawa potential energy of the pair (amu c^2).
  Double_t Potential(particle* p1_, particle* p2_);

  /// \brief Set the coupling g in amu c^2 fm.
  void     SetCoupling(Double_t g_){g = g_;}
  /// \brief The coupling g in amu c^2 fm.
  Double_t GetCoupling(){return g;}
  /// \brief Set the range a in fm.
  void     SetRange(Double_t a_){a = a_;}
  /// \brief The range a in fm.
  Double_t GetRange(){return a;}
protected:
  Double_t g = 0.1;    ///< coupling in amu c^2 fm
  Double_t a = 1.4138; ///< range in fm
private:
};

/// \brief Woods-Saxon potential (nucleus-nucleus optical potential, real part).
///
/// \f[ U = -\frac{V_0}{1 + e^{(r-R)/a}} \f]
///  - V0 [amu c^2] : depth. V0 > 0 is attractive. Default 0.0536772 amu c^2 = 50 MeV.
///  - R [fm] : radius. Default 8.886 fm = \f$1.2\,(4^{1/3} + 197^{1/3})\f$ fm, alpha + Au.
///  - a [fm] : diffuseness. Default 0.65 fm.
///
/// The force is finite everywhere; at r = 0 its direction is undefined and the null vector
/// is returned. It acts on every pair of particles of the event (independent of charge and
/// mass), so with many bodies each pair feels the nuclear potential; the same R is used for
/// every pair.
class woodsSaxonForce: public force{
public:
  /// \brief Make the force.
  /// \param V0_ depth in amu c^2 (default 0.0536772 = 50 MeV)
  /// \param R_ radius in fm (default 8.886)
  /// \param a_ diffuseness in fm (default 0.65)
  woodsSaxonForce(Double_t V0_ = 0.0536772, Double_t R_ = 8.886, Double_t a_ = 0.65):force(){
    V0 = V0_; R = R_; a = a_;
  }
  /// \brief Woods-Saxon force on p1_ by p2_ (amu c^2/fm).
  TVector3 Evaluate(particle* p1_, particle* p2_);
  /// \brief Woods-Saxon potential energy of the pair (amu c^2).
  Double_t Potential(particle* p1_, particle* p2_);

  /// \brief Set the depth V0 in amu c^2.
  void     SetDepth(Double_t V0_){V0 = V0_;}
  /// \brief The depth V0 in amu c^2.
  Double_t GetDepth(){return V0;}
  /// \brief Set the radius R in fm.
  void     SetRadius(Double_t R_){R = R_;}
  /// \brief The radius R in fm.
  Double_t GetRadius(){return R;}
  /// \brief Set the diffuseness a in fm.
  void     SetDiffuseness(Double_t a_){a = a_;}
  /// \brief The diffuseness a in fm.
  Double_t GetDiffuseness(){return a;}
protected:
  Double_t V0 = 0.0536772; ///< depth in amu c^2
  Double_t R  = 8.886;     ///< radius in fm
  Double_t a  = 0.65;      ///< diffuseness in fm
private:
};

#endif
