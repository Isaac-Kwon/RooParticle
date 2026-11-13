/// \file experiment.hpp
/// \brief Event factory that stamps out events from a template (\ref experiment, \ref experimentP).

#ifndef H_EXPERIMENT
#define H_EXPERIMENT 1;

#include "TROOT.h"
#include "particle.hpp"
#include "event.hpp"

#include "iostream"
#include "cstdarg"

/// \brief Makes a fresh event for each run from a template event.
///
/// Template pattern: configure one event once (integrator, forces, volume,
/// time step, ...) and give it with SetTemplate(). MakeEvent() then copies
/// that template (the event copy constructor) and adds the particles of the
/// run to the copy. The experiment keeps only the latest event, reachable
/// with GetEvent().
///
/// Ownership: the template is not owned (the caller keeps it alive while the
/// experiment is used). The created event is heap-allocated by MakeEvent()
/// and the experiment holds the pointer but never deletes it by itself: call
/// DelEvent() when the event is no longer needed (it also deletes the
/// particles and inspectors added to the event, see event). The particles
/// passed to MakeEvent() are owned by the event from then on, so create them
/// with \c new and do not pass one particle to two events. Calling MakeEvent()
/// again without DelEvent() leaks the previous event.
///
/// \warning The default constructor (also used by experimentP) leaves the template and event pointers
/// uninitialised; use experiment(event*) or call SetTemplate() before
/// MakeEvent(). GetEvent() returns a null pointer only after the
/// experiment(event*) constructor (and before the first MakeEvent()).
/// DelEvent() does not reset the pointer: GetEvent() afterwards returns a
/// dangling pointer, and a second DelEvent() deletes it twice.
class experiment{
public:
  /// \brief Empty experiment; the template must be set with SetTemplate().
  experiment(){}
  /// \brief Experiment with a template.
  /// \param eventT_ the template event (not owned, must outlive the experiment)
  experiment(event* eventT_);
  /// \brief Set (or replace) the template event.
  /// \param eventT_ the template event (not owned)
  void SetTemplate(event* eventT_){eventT=eventT_;}
  /// \brief Create a new event as a copy of the template and add particles to it.
  ///
  /// Variadic C call: \c ninparticle_ must equal the number of \c particle*
  /// arguments that follow. There is no type or count check, a wrong count
  /// reads garbage from the stack (undefined behaviour). Each argument must
  /// be exactly of type \c particle* (cast a pointer to a derived class
  /// first); passing \c 0 / \c NULL literals as an int is also wrong.
  /// The new event replaces the previous one (the previous is not deleted).
  /// \param ninparticle_ number of particle pointers that follow
  /// \param ... the particles, added in this order; owned by the new event
  void MakeEvent(Int_t ninparticle_, ...);
  /// \brief Mark the current event as off. The event itself is kept; only an internal flag changes (nothing else reads it).
  void OffEvent(){eventOn=false;}
  /// \brief Mark the current event as off and delete it (with its particles and inspectors).
  /// \warning The stored pointer is not reset, see the class description.
  void DelEvent(){eventOn=false; delete currentEvent; currentEvent=nullptr;}
  /// \brief The event created by the last MakeEvent() (owned by the experiment's user: see DelEvent()).
  event* GetEvent(){return currentEvent;}

  friend class experimentP;
protected:
private:
  event * eventT = nullptr;    ///< template event, not owned
  Bool_t eventOn = false;      ///< true after MakeEvent(), false after OffEvent()/DelEvent(); informational only
  event * currentEvent = nullptr; ///< last created event (heap)
};

/// \brief experiment with a fixed particle count: MakeEvent() takes no count.
///
/// Same template pattern and ownership as \ref experiment. The number of
/// particles per event is stored with the object (default 1), so the call is
/// \c MakeEvent(p1, p2, ...) with exactly that many pointers.
class experimentP: public experiment{
public:
  /// \brief Experiment that builds events of \c ninparticle_ particles.
  /// \param ninparticle_ number of particles per event
  experimentP(Int_t ninparticle_){ninparticle=ninparticle_;}
  /// \brief Change the number of particles per event.
  /// \param ninparticle_ number of particles per event
  void SetNInParticle(Int_t ninparticle_){ninparticle = ninparticle_;}
  /// \brief Number of particles per event.
  Int_t GetNInParticle(){return ninparticle;}
  /// \brief Create a new event from the template and add exactly GetNInParticle() particles.
  ///
  /// The first particle is always added (even if the count is 0 or negative);
  /// the remaining count-1 are read from the variadic list. The number of
  /// pointers passed must equal GetNInParticle(): there is no check, a
  /// shortfall is undefined behaviour. This hides experiment::MakeEvent.
  /// The new event replaces the previous one (not deleted); the particles are
  /// owned by the new event.
  /// \param p1 first particle (must be non-null)
  /// \param ... the remaining particle* arguments
  void MakeEvent(particle* p1, ...);
private:
  Int_t ninparticle=1; ///< particles per event
};


#endif