/// \file RPConfig.hpp
/// \brief Compile-time limits and default tolerances (macros).

/// Number of steps that event::DeriveMAX always takes (it is a fixed count, not a safety
/// cap: DeriveMAX runs exactly this many steps). Unit: steps. Used only by DeriveMAX in
/// src/event.cpp. Unrelated to the stop reason named "MAXSTEPS" (event::kStopMaxSteps),
/// which concerns MAXINSPECTSTEPS.
#define MAXSTEPS 100000
/// Safety cap of event::DeriveInspect (number of steps). 0 means unlimited.
/// Unit: steps. Initial value of event's \c maxInspectSteps; the test programs also
/// compare their own step counter with it.
#define MAXINSPECTSTEPS 100000000
/// Default relative tolerance of the adaptive integrator event::kRK45
/// (event::SetTolerance): relative to the size of each position / velocity.
/// Dimensionless. Initial value of event's \c rk45Rel.
#define RK45MAXSUBSTEPS 1000000 ///< Substeps event::kRK45 may try in one DeriveDT() call before it gives up (event::HasIntegratorFailed()).
#define RK45RELTOL 1e-8
/// Default absolute tolerance of the adaptive integrator event::kRK45
/// (event::SetTolerance): fm for positions, c for velocities.
/// Initial value of event's \c rk45Abs.
#define RK45ABSTOL 1e-10
