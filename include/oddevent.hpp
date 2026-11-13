#ifndef H_ODDEVENT
#define H_ODDEVENT 1

#include "TROOT.h"
#include "event.hpp"

/// \brief Variant of event which moves each particle right after its force is calculated.
///
/// The particles are visited in the order of addition: the force on a particle is evaluated
/// from the present positions of all the others (every ordered pair, every force), the
/// particle is moved at once, and the later particles therefore see the already moved ones.
/// An invincible particle is not moved but counts the step.
///
/// It has its own scheme (TAYLOR-like, non-relativistic) and ignores the integrator,
/// event::SetPairSymmetric() and event::SetRelativistic() settings.
///
/// Many bodies: because the update is sequential, action and reaction are not applied at
/// the same state. As soon as more than one particle moves (two or more moving bodies, the
/// usual case in a many-body event), the total momentum is not exactly conserved, and the
/// energy error is larger than with the simultaneous schemes of event. With a single moving
/// particle (all others invincible) it is equivalent to the TAYLOR scheme. A copy of an
/// oddevent made with `new event(*ev)` is a plain event, not an oddevent.
class oddevent : public event{
    public:
    /// \brief Advance every particle by dt, one after another.
    /// \param dt time step in fm/c
    ///
    /// Each non-invincible particle receives the sum of the forces from all other particles
    /// (positions as left by the particles moved before it) and moves with
    /// `x += v dt + a dt^2/2 ; v += a dt`. Every particle counts one step and moving ones
    /// advance their time by dt.
    void DeriveDT(Double_t dt = 1);
    protected:
    /// \brief Always kFALSE: this scheme never uses the relativistic option.
    Bool_t UseRelativistic(){return kFALSE;}
};

#endif
