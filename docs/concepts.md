\page concepts Concepts

This page explains the model behind the \ref tutorial "Tutorial": what the
objects are and who deletes them, what one step does, what each integrator
computes and conserves, exactly when an event stops, why a finite start and
stop distance matters and what to do about it, and what makes a run
reproducible. Everything here was checked against the source code; where the
root `README.md` or `TODO.md` says otherwise, the code is described.

The code is built for events of **three or more bodies** interacting
pairwise. Two-body events are the special case with an exact solution, used
for validation (\ref validation "Validation").

# Units and constants {#con_units}

| Quantity | Unit | Note |
|---|---|---|
| length | fm | \f$ 10^{-15} \f$ m |
| mass | amu | rest mass, `particle::GetM0()` |
| velocity | c | speeds are fractions of c, \f$ c = 1 \f$ |
| charge | e | elementary charge |
| time | fm/c | \f$ 3.33564 \times 10^{-24} \f$ s |
| energy | amu c\f$^2\f$ | 1 amu c\f$^2\f$ = 931.494 MeV |
| momentum | amu c | 931.494 MeV/c |
| force | amu c\f$^2\f$/fm | 931.494 MeV/fm |

Conversions you will need:

- kinetic energy in MeV: \f$ T = \tfrac12 m v^2 \times 931.494 \f$ (with m in
  amu, v in c); a 4.66 MeV alpha has \f$ v = 0.05\,c \f$.
- speed from a kinetic energy: \f$ v = \sqrt{2T/(m c^2)} \f$
  non-relativistically, \f$ v = \sqrt{1 - (1 + T/(m c^2))^{-2}} \f$ exactly.

**The Coulomb constant.** \f$ U = k q_1 q_2 / r \f$ with
\f$ k = e^2/(4\pi\varepsilon_0) = \hbar c\,\alpha = 1.439965 \f$ MeV fm.
Divided by 931.494 MeV per amu c\f$^2\f$ this is
\f$ k = 0.00154586 \f$ amu c\f$^2\f$ fm/e\f$^2\f$. The code uses
`coulombK = 0.00154629` (`coulombForce::GetK()`), which is 1.44036 MeV fm with
931.494 MeV/amu, i.e. \f$ 2.8 \times 10^{-4} \f$ larger (it corresponds to
1.43996 MeV fm only with 931.24 MeV/amu). Every Coulomb result therefore
carries this relative bias in \f$ \kappa \f$; for a fixed `b` and `v` the
Rutherford angle changes by about \f$ 2.8\times10^{-4}\sin\theta \f$ rad.
`SetK()` changes the value per force object.

# Object model and ownership {#con_ownership}

| Class | What it is |
|---|---|
| \ref particle | point mass: position, velocity, rest mass, charge (0 for a plain particle), step count, evolved time, optional track |
| \ref EMparticle | a particle constructed with a charge; adds nothing else |
| \ref force | pairwise interaction: `Evaluate(p1, p2)` (force on p1 by p2) and `Potential(p1, p2)` |
| \ref event | one collision: particles, forces, inspectors, integrator settings |
| \ref oddevent | event with a sequential update (section \ref con_odd) |
| \ref inspector, \ref inspectorP, \ref inspectorE | terminating conditions (section \ref con_stop) |
| \ref eventVolume, \ref eventVGeneral | box used by the `VOUT`/`VOUTA` inspectors |
| \ref path | track of a recording particle |
| \ref asymptote | exact two-body asymptotes from one state (section \ref con_asymptote) |
| \ref asymptoteN | asymptotic velocities of any number of bodies, with error estimates (section \ref con_nbody_asymptote) |
| `runinfo` | run metadata written into the ROOT file |

An **invincible** particle (constructor flag) is never moved: it gets no
force, its velocity is not changed, and its time stays 0, but it counts steps
like the others and acts on every other particle. It models an infinitely
heavy fixed centre.

**Ownership.** The event deletes, in its destructor:

- the particles given to `AddParticle()`,
- the inspectors given to `AddInspector()` (also ones added after a copy),
- its pre-derived copy (`preDerive()`).

It does **not** delete the forces (`AddForce()`) nor the volume
(`SetVolume()`): they are shared between an event and its copies and may be
objects on the stack. Consequences: particles and inspectors must be created
with `new`; never give one particle or inspector to two events; keep each
force and volume alive while any event that uses it exists. Assignment
(`operator=`) is deleted.

## What copying an event copies {#con_copy}

`event(event&)`, the copy constructor used for the template-and-copies
pattern:

- **particles**: each is copied with `new particle(*p)`: position,
  velocity, mass, charge, invincible flag, **step count and evolved time**,
  record flag, record period and the recorded track. The copy is always a
  plain `particle`: a subclass (including `EMparticle`, which has no data of
  its own, or your own subclass with extra members or virtual functions)
  loses its type. The charge survives because it is a member of `particle`.
- **forces**: the same pointers (shared, not copied).
- **inspectors**: new inspectors with the same method, value and abort flag.
  Pair inspectors are re-linked to the copied particles with the same index;
  event inspectors watch the copy. An `EDG`/`EDL` reference energy is not
  copied (taken again at the first inspection).
- **settings**: integrator, inspect mode, step cap, pair symmetry,
  relativistic flag, RK45 tolerances, volume pointer.
- **not copied**: the RK45 substep suggestion and its counters, the last stop
  reason, the pre-derived copy, the Verlet force cache.

Because the step count is copied, `CNT` inspectors of a copy of a derived
event continue from the source's count. A copy of an `oddevent` made with
`new event(*odd)` is an `event` (with the integrator setting of the source,
`TAYLOR` by default); use `new oddevent(*odd)`.

If a pair inspector of the source refers to a particle that is not in the
source event, it is not copied (with a message on stderr): the copy then has
one terminating condition less. Add inspectors only for particles of the same
event.

# One step {#con_step}

`event::DeriveDT(dt)` advances every non-invincible particle by `dt` with the
selected integrator, then closes the step for **every** particle: step count
+1, evolved time + `dt` (not for invincible ones), and a track record if the
particle records and its step count is a multiple of its record period.
`DeriveDTN(n, dt)` is `n` such calls; `DeriveMAX(dt)` makes 100,000 of them;
`DeriveInspect()` repeats them until the inspectors stop it (section
\ref con_stop).

Forces depend on positions only (and on particle constants such as charge and
rest mass): there is no time dependence and no velocity dependence. The
integrators except `TAYLOR` and `oddevent` evaluate all forces at one common
state before moving any particle.

## Pair-symmetric force evaluation {#con_pair}

With `SetPairSymmetric(kTRUE)` (default), each unordered pair is evaluated
once, `Evaluate(p1, p2)` goes to p1 and its negative to p2 (Newton's third
law). A pair of two invincible particles is skipped; a pair with one
invincible particle is evaluated once, for the moving one. With `kFALSE`,
every ordered pair is evaluated separately: use that for a force whose
`Evaluate(p2, p1)` is not `-Evaluate(p1, p2)`. The setting is used by
`SEULER`, `VERLET`, `RK4`, `RK45` and by `GetMaxAcceleration()` (`ACG`/`ACL`).
`TAYLOR` and `oddevent` always evaluate every ordered pair.

For a correctly written force the two settings give the same trajectory to
round-off (`test_event_check`: \f$ 10^{-8} \f$ fm over 10000 steps) with 9/12
of the evaluations for 3 moving and 2 fixed particles. Pair symmetry makes
the total momentum conserved to round-off in every integrator.

# Integrators {#con_integrators}

Selected with `SetIntegrator("NAME")` or `SetIntegrator(event::kName)`. Below,
\f$ a(x) = F(x)/m_0 \f$ and the step is `dt` = \f$ h \f$.

**TAYLOR** (`kTaylor`, **default**, the original scheme):
\f[ x_{n+1} = x_n + v_n h + \tfrac12 a(x_n) h^2, \qquad v_{n+1} = v_n + a(x_n) h. \f]
First order (the velocity update is Euler's), not symplectic: the energy
drifts linearly. One evaluation per ordered pair with a moving first
particle. Ignores `SetRelativistic()` and `SetPairSymmetric()`.

**SEULER** (`kSymplecticEuler`):
\f[ v_{n+1} = v_n + a(x_n) h, \qquad x_{n+1} = x_n + v_{n+1} h. \f]
First order, symplectic: the energy error stays bounded, angular momentum is
conserved exactly for central forces. One force evaluation per step.

**VERLET** (`kVelocityVerlet`):
\f[ x_{n+1} = x_n + v_n h + \tfrac12 a(x_n) h^2, \qquad v_{n+1} = v_n + \tfrac12\left(a(x_n) + a(x_{n+1})\right) h. \f]
Second order, symplectic and time-reversible; exact angular momentum for
central forces. The force at the end of a step is reused at the start of the
next one, so it costs **one** evaluation per step; the reuse is skipped
(two evaluations) if any particle was moved from outside (`applyX()` etc.)
since the last step. The cache does **not** notice a change of a force
parameter (`SetK()`, `SetSoftening()`, ...) between steps: the first step
after such a change starts from the old force. Change parameters only before
deriving, or call `SetIntegrator()` again (which clears the cache). In the relativistic mode it is the kick-drift-kick form
on the momentum.

**RK4** (`kRK4`): classical fourth-order Runge-Kutta on the state
\f$ (x, u) \f$ with \f$ u = v \f$ (or \f$ u = p \f$ relativistic),
\f$ \dot x = v,\ \dot u = F/m_0 \f$ (or \f$ F \f$). Four evaluations per
step. Not symplectic: energy and angular momentum drift slowly.

**RK45** (`kRK45`), adaptive, see \ref con_rk45.

Measured by `test_integrator_compare` (alpha on fixed Au, b = 50 fm, from
1000 fm, fixed time 40000 fm/c, `dt` = 8 fm/c, reference RK4 at
`dt` = 0.125): relative energy change and error of the final direction:

| Name | Order | Evaluations per step | Energy change | Direction error (deg) |
|---|---|---|---|---|
| `TAYLOR` | 1 | 1 per ordered pair | 8.7e-4 | 3.2e-2 |
| `SEULER` | 1 | 1 | 2.0e-5 | 2.3e-4 |
| `VERLET` | 2 | 1 | 2.9e-10 | 3.2e-5 |
| `RK4` | 4 | 4 | 1.2e-15 | 7.0e-11 |

The direction error at a fixed time includes the error of the position along
the orbit. The error of the orbit itself (change of the asymptotic angle
during a run) is much smaller for VERLET: \f$ 1.5\times10^{-12} \f$ degree at
`dt` = 1 (\ref tut_accuracy). Measured conservation for two free charges over
60000 fm/c (`test_event_check`): momentum to \f$ 10^{-14} \f$ for all event
integrators; angular momentum \f$ 2\times10^{-14} \f$ (VERLET, SEULER) and
\f$ 4\times10^{-14} \f$ (RK4, RK45 over that run, not exact in general);
`oddevent` \f$ 1.5\times10^{-6} \f$ in momentum.

Which to use: `VERLET` for fixed steps and long runs; `RK45` when the
distances vary over orders of magnitude (approach, encounter, flight to large
distance), which is the normal case for scattering; `RK4` for a short smooth
run at a fixed step. `TAYLOR` stays the default only so that old results are
reproduced (`TODO.md` lists the decision as open).

## RK45 in detail {#con_rk45}

Dormand-Prince 5(4) with the first-same-as-last property: the fifth-order
solution is kept, the embedded fourth-order one gives the error estimate.
`DeriveDT(dt)` advances exactly `dt` and counts one step, taking substeps
\f$ h \le dt \f$ inside. A substep is accepted when the RMS over the moving
particles of the scaled errors
\f[ \frac{|\delta x|}{\mathrm{abs} + \mathrm{rel}\,|x|}, \qquad \frac{|\delta u|}{\mathrm{abs} + \mathrm{rel}\,|u|} \f]
is at most 1 (positions in fm, velocities in c; momentum divided by the rest
mass in the relativistic mode). Defaults \f$ \mathrm{rel} = 10^{-8} \f$,
\f$ \mathrm{abs} = 10^{-10} \f$ (`RPConfig.hpp`); `SetTolerance(rel, abs)`
refuses zero or negative values. The next substep is
\f$ h \cdot \min(5, \max(0.2,\ 0.9\,\mathrm{err}^{-1/5})) \f$, not enlarged
right after a rejection, and the last suggestion is reused by the next call.

Cost: 6 evaluations per substep plus 1 per call. Practical consequences:

- The call `dt` caps the substep. With a `dt` smaller than the substep the
  tolerance would allow, every call is one substep and the tolerance has no
  effect. Under `DeriveInspect()` the `dt` is also the inspection interval.
- The relative tolerance is relative to \f$ |x| \f$, the distance from the
  **coordinate origin**, not to the pair separation. Put the origin at the
  collision.
- Below \f$ 10^{-12}\,dt \f$ a substep is accepted anyway with a warning
  (once per event), whatever its error.
- A singular force (unsoftened attraction to \f$ r = 0 \f$) makes the
  substeps shrink inside one call, where inspectors and the step cap are not
  looked at. After `RK45MAXSUBSTEPS` (\f$ 10^6 \f$, `RPConfig.hpp`) substeps
  in one call the integrator gives up: it restores the last accepted state,
  `HasIntegratorFailed()` becomes true and stays true, and `DeriveInspect()`
  ends with `kStopFailed`. The "last accepted state" may itself be broken (in
  the softening example of the tutorial it is NaN, because substeps at the
  minimum size are accepted whatever their error): discard such an event.

# oddevent {#con_odd}

`oddevent::DeriveDT()` loops over the particles; for each moving particle it
sums the forces of all others at their **current** positions and moves it at
once with the `TAYLOR` update. Later particles see earlier ones already
moved. With one moving particle this is identical to `event` with `TAYLOR`.
With several, action and reaction are evaluated at different states: the
total momentum is not conserved (\f$ \sim 10^{-6} \f$ relative in the tests)
and the result depends on the order of the particles. It ignores the
integrator, pair-symmetry and relativistic settings.

# Relativistic mode {#con_rel}

With `SetRelativistic(kTRUE)` the integrated state is \f$ (x, p) \f$ with
\f[ p = \gamma m_0 v, \qquad \frac{dp}{dt} = F(x), \qquad v = \frac{p}{\sqrt{m_0^2 + p^2}} \quad (c = 1). \f]
It is used by `SEULER`, `VERLET`, `RK4` and `RK45`. `TAYLOR` and `oddevent`
silently stay non-relativistic.

What is integrated: relativistic kinematics of each particle in the
instantaneous (static) pairwise forces. What is ignored: retardation (the
force uses the other particle's present position), magnetic forces between
moving charges, radiation. For a charge in the field of a fixed centre this
is exact. For two moving charges the neglected terms are of order
\f$ (v/c)^2 \f$, the same order as the relativistic correction kept, so the
mode does not make such an event more accurate in a controlled way. Speeds
must stay below 1.

Energies: `GetNetKERel()` is \f$ \sum (\gamma - 1) m_0 \f$ and
`GetNetERel()` adds the potential energy; that sum is conserved in this mode.
`GetNetKE(kTRUE)` and `particle::GetKE(kTRUE)` are
\f$ \tfrac12 \gamma m_0 v^2 \f$, which is not the relativistic kinetic energy
and is not conserved. `GetNetMomentum()` and `GetNetAngularMomentum()` use
\f$ \gamma m_0 v \f$ when the mode is in effect. The inspectors `KEG`/`KEL`,
`PKG`/`PKL` and `EDG`/`EDL` use the relativistic kinetic energy
\f$ \sum(\gamma - 1) m_0 \f$ when the flag is set and the integrator is not
`TAYLOR` (they test the flag and the integrator, not the event class, so an
`oddevent` with the flag set and a non-`TAYLOR` integrator setting is
measured relativistically although it moves non-relativistically).
\ref asymptote is non-relativistic; \ref asymptoteN applies only its
first-order correction in this mode (section \ref con_nbody_asymptote).

# Terminating conditions {#con_stop}

## How DeriveInspect works {#con_stop_algo}

`DeriveInspect(iperiod = 1, dt = 1)`:

1. If any inspector has an undefined method: message, stop reason
   `kStopInvalid`, return `kFALSE`, nothing derived.
2. Loop with a counter `i` = 0, 1, 2, ... (steps of this call):
   if `i` is a multiple of `iperiod`, call `Inspect()`; if it says stop,
   end. Otherwise, if a cap is set and `i` has reached it, end with
   `kStopMaxSteps`. Otherwise `DeriveDT(dt)`; if the integrator failed in
   it (`HasIntegratorFailed()`), end with `kStopFailed`.
3. Return `kTRUE` only for `kStopInspector`.

So the inspectors are evaluated before the first step, at most
`SetMaxInspectSteps(n)` steps are made per call (default
100,000,000 = `MAXINSPECTSTEPS`, 0 = unlimited), and with `iperiod > 1` the
event stops at a multiple of `iperiod` steps of the call.

`Inspect()`:

1. Every **abort** inspector (`SetAbort()`) is evaluated in order; the first
   one satisfied stops the event with `kStopAbort`, whatever the others say.
2. The other inspectors are all evaluated and combined with **AND** (default,
   `kInspectAND`: all satisfied) or **OR** (`SetInspectMode(kInspectOR)`:
   any satisfied).
3. No normal inspectors: with abort inspectors, continue (only an abort or
   the cap ends the event, and `DeriveInspect()` returns `kFALSE`); with no
   inspectors at all, stop at once.

| `GetStopReason()` | `GetStopReasonName()` | Meaning | `DeriveInspect()` returns |
|---|---|---|---|
| `kStopNone` | `NONE` | not derived by `DeriveInspect()` yet | |
| `kStopInspector` | `INSPECTOR` | normal inspectors satisfied | `kTRUE` |
| `kStopAbort` | `ABORT` | an abort inspector fired | `kFALSE` |
| `kStopMaxSteps` | `MAXSTEPS` | the step cap was reached | `kFALSE` |
| `kStopInvalid` | `INVALID` | undefined method, nothing derived | `kFALSE` |
| `kStopFailed` | `FAILED` | the integrator gave up (RK45 on a singular force); the state is not usable | `kFALSE` |

All comparisons are strict (`>` or `<`). The method can be given as a string
or as the enumerator `inspector::kDEG` etc.; names are case-sensitive.

## Pair inspector methods (inspectorP) {#con_stop_pair}

\f$ \Delta x = x_1 - x_2 \f$, \f$ \Delta v = v_1 - v_2 \f$ of the two
particles given to the constructor.

| Method | Evaluated quantity | Satisfied when |
|---|---|---|
| `COS` | cosine of the angle between \f$ \Delta x \f$ and \f$ \Delta v \f$, clamped to [-1, 1] | value < criterion |
| `DEG` | that angle in degree | value < criterion |
| `RAD` | that angle in radian | value < criterion |
| `CNT` | the smaller step count of the two particles | value > criterion |
| `SPDG` | relative speed \f$ \lvert\Delta v\rvert \f$ | value > criterion |
| `SPDL` | relative speed | value < criterion |
| `DSTG` | distance \f$ \lvert\Delta x\rvert \f$ | value > criterion |
| `DSTL` | distance | value < criterion |
| `VRG` | radial velocity \f$ \Delta x\cdot\Delta v/\lvert\Delta x\rvert \f$ | value > criterion (0: receding) |
| `VRL` | radial velocity | value < criterion (0: approaching) |

The angle is \f$ \mathrm{atan2}(\lvert\Delta x \times \Delta v\rvert, \Delta x \cdot \Delta v) \f$,
accurate also near 0 and 180 degree where \f$ \arccos \f$ of the cosine is
not. **Zero-vector rule:** if \f$ \Delta x \f$ or \f$ \Delta v \f$ is exactly
the zero vector, the angle (and cosine) is undefined: `Evaluate()` returns NaN
and `COS`/`DEG`/`RAD` are not satisfied; if \f$ \Delta x = 0 \f$, the same
for `VRG`/`VRL`. Distance and speed methods work as usual.

## Event inspector methods (inspectorE) {#con_stop_event}

| Method | Evaluated quantity | Satisfied when |
|---|---|---|
| `KEG` / `KEL` | net kinetic energy of the moving particles (relativistic when the mode is in effect) | > / < criterion |
| `PEG` / `PEL` | net potential energy of all pairs | > / < criterion |
| `PKG` / `PKL` | absolute value of PE / KE, KE as for `KEG` (largest double if KE = 0) | > / < criterion |
| `EDG` / `EDL` | relative energy change abs((E - E0)/E0), E0 at the first evaluation (absolute change if E0 = 0) | > / < criterion |
| `TMG` / `TML` | evolved time, largest over the particles | > / < criterion |
| `DXG` | largest step count | > criterion |
| `DXL` | smallest step count | < criterion |
| `DNG` | largest step count | > criterion |
| `DNL` | smallest step count | < criterion |
| `DSG` / `DSL` | smallest distance between any two particles | > / < criterion |
| `ACG` / `ACL` | largest acceleration F/m0 of the moving particles | > / < criterion |
| `VOUT` | some moving particle is outside the volume | (criterion unused) |
| `VOUTA` | all moving particles are outside the volume | (criterion unused) |
| `ASYD` | largest remaining deflection of a moving particle (or cluster), rad, first order (`asymptoteN::RemainingDeflection()`) | < criterion **and** every interacting pair receding |
| `ASYV` | largest remaining relative velocity change \f$ \lvert v_\infty - v\rvert / \lvert v_\infty\rvert \f$ (`asymptoteN::RemainingChange()`) | < criterion **and** every interacting pair receding |

Notes: `DXG` and `DNG` test the same thing, as do `DXL` and `DNL`, but
`Evaluate()` of `DX*` returns the largest and of `DN*` the smallest count,
so the printed value of `DXL` and `DNG` is not the one compared. Without a
volume (`SetVolume()`), every particle counts as inside: `VOUT`/`VOUTA` are
never satisfied. `DSG`/`DSL` and the potential energy include pairs of two
invincible particles. `ResetReference()` forgets E0. `ASYD`/`ASYV` are
relative to each particle's own speed, so a slow particle (a recoiling heavy
target) dominates them; `Evaluate()` returns the remaining quantity without
the receding test. They compute a first-order (`kBorn`) solution at every
inspection: a few microseconds for the plain Coulomb and gravity forces,
much more for any other force class (numerical quadrature), including a
subclass of `coulombForce`.

## Which inspectors fit which event {#con_stop_many}

- One projectile on a fixed centre: `DEG < 1 AND CNT > N` (the project's
  default; the step minimum keeps near head-on events from stopping at the
  turning point), or `VRG > 0 AND DSTG > R`. With \ref asymptote the stop
  can be anywhere after the collision.
- Several moving bodies: `ASYD` (or `ASYV`) at a value \f$ \epsilon \f$,
  followed by \ref asymptoteN. The error of the pairwise result then falls
  like \f$ \epsilon^2 \f$ (in the tutorial example \f$ 0.04\,\epsilon^2 \f$
  relative to the speed; \ref con_nbody_asymptote).
  Alternatives that say nothing about the remaining change: one `VRG > 0` per
  pair with `PKL` or `DSG`. `DEG` on one pair is not a criterion for the
  event.
- Always: an abort on `EDG` (numerical failure) and, for attractive forces,
  on `DSL` (close encounter), plus a step cap suited to the event.

Stopping "far enough" never means "at infinity": see the next section.

# The finite-distance problem and the two-body asymptote {#con_asymptote}

A Coulomb interaction never ends: at distance \f$ r \f$ the potential energy
of a pair is \f$ \kappa/r \f$, the remaining deflection of the order of
\f$ \kappa/(\mu v^2 r) \f$. The simulation must start and stop at finite
distance, so the state at the start is not the incoming asymptote and the
state at the stop is not the outgoing one. For alpha on Au at 0.05 c from
1000 fm: \f$ v_\infty = 0.0512 \f$ c instead of 0.05 c, effective
\f$ b = 48.82 \f$ fm instead of 50 fm, and the Rutherford angle from the
nominal values is 1.07 degree too large (\ref tut_physics).

**Two bodies, pure 1/r.** For \f$ U = \kappa/r \f$ the relative orbit is a
conic section fixed by three conserved quantities of the relative coordinates
\f$ r = x_1 - x_2 \f$, \f$ v = v_1 - v_2 \f$, with reduced mass
\f$ \mu = m_1 m_2/(m_1+m_2) \f$ (or \f$ \mu \f$ = the moving mass if the other
body is invincible):
\f[ E = \tfrac12 \mu v^2 + \frac{\kappa}{r}, \qquad L = \mu\, r \times v, \qquad A = (\mu v) \times L + \mu \kappa\, \hat r . \f]
\f$ A \f$ (Laplace-Runge-Lenz) points from the centre to the closest
approach. From any one state, \ref asymptote computes
\f[ v_\infty = \sqrt{2E/\mu}, \quad e = \sqrt{1 + \frac{2 E L^2}{\mu \kappa^2}}, \quad \sin\frac{\theta}{2} = \frac1e, \quad b = \frac{\lvert L\rvert}{\mu v_\infty}, \f]
the incoming and outgoing relative velocities, the closest approach
\f$ r_{min} = a(1+e) \f$ with \f$ a = \lvert\kappa\rvert/(\mu v_\infty^2) \f$
(repulsive), the centre-of-mass velocity and the lab velocities and angles of
both bodies. The result satisfies the Rutherford relation
\f$ \tan(\theta/2) = \kappa/(\mu v_\infty^2 b) \f$ exactly
(`asymptote::RutherfordAngle()`); `test_asymptote_check` shows it is
invariant along a trajectory to \f$ 10^{-13} \f$.

`GetStatus()`: `kUnbound` (hyperbola), `kHeadOn` (\f$ L = 0 \f$, repulsive,
\f$ \theta = 180 \f$ degree, \f$ r_{min} = \kappa/E \f$), `kFree`
(\f$ \kappa = 0 \f$), `kRadialAttractive`, `kParabolic`, `kBound` (no
asymptote), `kInvalid` (zero distance, non-positive mass, both bodies fixed).
Only the first three have asymptotes (`HasAsymptote()`); otherwise the
asymptotic getters return 0, never NaN. `kRadialAttractive` is returned for
any radial motion under an attraction, also for a pair that recedes with
positive energy and therefore does have an outgoing asymptote;
\ref asymptoteN handles that case.

**Validity.** The result is exact only if all of these hold: exactly two
bodies; a pure \f$ 1/r \f$ interaction (no nuclear force, screening or
softening); non-relativistic mechanics; and the given state itself is
correct. A short-range extra force is allowed if the state is taken where it
is negligible. Otherwise the result is an approximation that improves with
the distance at which the state is taken.

**Lab angles of a free target.** "At rest" means at rest on the asymptote. A
target placed at rest at a finite distance from a charged projectile is
already being pushed, so on its incoming asymptote it moves; `GetThetaLab2()`
then measures from that motion. Shift all velocities by `GetV2In()` to make
the target at rest at infinity (\ref tut_forces).

## Three or more bodies {#con_nbody}

There is no closed form, and the Coulomb tail makes "integrate until it
stops changing" expensive: the remaining velocity change falls only like
\f$ 1/R \f$ with the distance \f$ R \f$. Two ways out:

- **Integrate far and extrapolate** (\ref tut_many, part B): integrate with
  `RK45` until every pair recedes and the closest pair is `R` apart, record
  the velocities, raise `R` tenfold, continue; extrapolate
  \f$ v_\infty \approx v(R) + (v(R) - v(R/10))/9 \f$. `asymptoteN::kBrute`
  does this efficiently (doubling time blocks) and is the reference.
- **Correct analytically from a moderate distance** with \ref asymptoteN
  (`kPairwise`), below. In the tutorial example this gives at 1414 fm the
  accuracy that the extrapolation reaches at \f$ 10^5 \f$ fm and the raw
  velocities at \f$ 10^7 \f$ fm (\ref tut_nbody_asymptote).

In both cases check energy, momentum (no invincible particles) and angular
momentum (about the fixed centre if there is one). Conservation laws show that
the integration is right, not that the stop was far enough.

## The N-body asymptote: asymptoteN {#con_nbody_asymptote}

\ref asymptoteN returns, for every moving particle of an event, the
asymptotic velocity \f$ v_{\infty} \f$ and an error estimate, from the present
state, for \f$ t \to +\infty \f$ (`kOutgoing`) or \f$ t \to -\infty \f$
(`kIncoming`, by time reversal: all velocities of a copy flipped). It is
about velocities only: with \f$ 1/r \f$ forces the positions never settle on
straight lines (the logarithmic Coulomb phase), the velocities converge like
\f$ 1/t \f$.

Notation for a pair \f$ a, b \f$: \f$ r = x_a - x_b \f$, \f$ w = v_a - v_b \f$,
\f$ R = \lvert r\rvert \f$, \f$ \mu_{ab} \f$ the reduced mass (the moving mass
if the partner is invincible), \f$ E_{ab} = \tfrac12 \mu_{ab} w^2 + U_{ab}(R) \f$.

**First order, `kBorn`.** Every particle moves on its present straight line;
the impulse of each pair force is integrated from now to infinity. For
\f$ U = \kappa/r \f$ (unsoftened Coulomb and gravity) in closed form,
\f[ I_{ab} = \int_0^\infty F_{ab}(r + w t)\,dt = \frac{\kappa}{w R}\left[\frac{\vec b}{R + s} + \hat w\right], \qquad s = r\cdot\hat w,\quad \vec b = r - s\,\hat w, \f]
for other forces by adaptive Gauss-Legendre quadrature of
`force::Evaluate()` along the line. Then
\f$ v_{\infty,a} = v_a + \sum_b I_{ab}/m_a \f$. The neglected curvature of
the paths makes the error \f$ O(\eta^2) \f$ per pair, with the strength
\f[ \eta_{ab} = \frac{\max(\lvert U_{ab}\rvert,\ R\,\lvert F_{r}\rvert)}{2 E_{ab}} \sim \frac{1}{R}, \f]
i.e. \f$ O(1/R^2) \f$. This is what `RemainingDeflection()`,
`RemainingChange()` and the inspectors `ASYD`/`ASYV` evaluate (the angle,
or the relative size, of \f$ v_\infty - v \f$).

**Pairwise, `kPairwise` (default).** Every pair is solved exactly as an
isolated two-body problem: with \ref asymptote for a pure \f$ 1/r \f$ pair
(closed form), otherwise by the orbit integral
\f[ \Delta\varphi = \int \frac{h\,dr}{r^2 \sqrt{2(E-U(r))/\mu - h^2/r^2}} \f]
(adaptive quadrature through the turning point, so an approaching pair,
also one passing through a nucleus, is solved too). This gives each pair's
relative velocity change \f$ \Delta w_{ab} = w_{ab,\infty} - w_{ab} \f$, and
the changes are added with the two-body momentum shares:
\f[ v_{\infty,a} = v_a + \sum_b \frac{\mu_{ab}}{m_a}\,\Delta w_{ab}. \f]
Exact for two bodies (\f$ 10^{-16} \f$ closed form, \f$ 10^{-15} \f$
quadrature) and momentum-conserving exactly. For three or more bodies the
missing piece is the **cross term**: pair \f$ ab \f$ is solved with its
present relative velocity, but the other bodies will still change it by
\f$ \Delta w^{ext}_{ab} = (\Delta v_a - \Delta v_b) - \Delta w_{ab} \sim 1/R \f$,
and the pair's own change responds to that at relative strength
\f$ \eta_{ab} \sim 1/R \f$: an error \f$ O(1/R^2) \f$, like `kBorn`, but with
a 10 to 60 times smaller coefficient in the tests of `asymptoteN.hpp`. The
energy at infinity (`GetEnergyInf()`) differs from the present energy
(`GetEnergyNow()`) by the same order.

**Error estimate** (`GetError(i)`, in c, a magnitude, not a strict bound):
\f[ \delta v_a = C \sum_b \frac{\mu_{ab}}{m_a}\, \eta_{ab}\, \lvert\Delta w^{ext}_{ab}\rvert \;+\; \text{quadrature error} \;+\; \text{cluster multipole error}, \qquad C = 2 \f]
(`SetCrossFactor()`; \f$ \eta \f$ is taken at the present distance, or at the
closest approach still ahead for an approaching pair). For two bodies only the
quadrature and round-off terms remain. Measured safety factor: actual /
estimate \f$ \le 0.23 \f$ for `kPairwise` and `kAuto` in the tests of
`asymptoteN.hpp` (3 and 4 bodies, fixed and free targets, Coulomb and
Coulomb + Woods-Saxon, stops from 15 to 23000 fm), 0.15 for the largest
error in the tutorial example: the estimate is 4 to 9 times too large. For `kBorn` the
estimate adds \f$ \lvert \Delta v_{Born} - \Delta v_{Pairwise}\rvert \f$ and is
tight (ratio \f$ \le 0.99 \f$).

**`kAuto` and `kBrute`.** `kAuto` continues a private copy with RK45 in
doubling time blocks until the pairwise estimate is below
`SetAutoTolerance()` (default \f$ 10^{-9} \f$ relative). `kBrute` continues a
copy in doubling blocks \f$ t_k \f$ and extrapolates with the \f$ 1/t \f$ law,
\f$ v_\infty = (t_2 v_2 - t_1 v_1)/(t_2 - t_1) \f$, until two extrapolations
agree within `SetBruteTolerance()` (\f$ 10^{-11} \f$); it never uses the
analytic methods and is the reference (a few hundred to a few thousand force
sweeps, to \f$ t \sim 10^9 \f$ fm/c). The caller's event is never changed.

**Clusters.** Bodies are merged while some pair has a negative pair energy in
its own centre-of-mass frame. A cluster is treated as one body (total mass and
charge at its centre of mass); every member gets the cluster's velocity (flag
`kInCluster`). Its internal motion and its multipole interaction with the
others (e.g. the dipole of a bound \f$ +q, -q \f$ pair) are not in
\f$ v_\infty \f$; the multipole part enters the error estimate.
`GetClusterIsolation()` compares the coupling to the other bodies with the
binding; above `SetIsolation()` (\f$ 10^{-2} \f$) the flag `kNotIsolated` says
that the binding may still be broken or made. A body bound to an invincible
particle has no asymptote (`kBoundToFixed`, \f$ v_\infty = 0 \f$).

**Flags.** `kApproaching` (a pair still approaches: handled, estimate less
certain), `kStrong` (some pair has \f$ C\eta > 0.5 \f$ and is disturbed by
others: not trustworthy), `kInCluster`, `kNotIsolated`, `kBoundToFixed`,
`kRelativistic`, `kNoPairwise` (non-central force or
`SetPairSymmetric(kFALSE)`: first order only), `kMarginal` (pair energy near
zero), `kSingular` (radial fall into an attractive centre), `kNotConverged`
(`kAuto`/`kBrute` hit `SetMaxBlocks()`), `kInvalid`. `IsReliable()` is false
when any of `kStrong`, `kNotIsolated`, `kMarginal`, `kSingular`,
`kNotConverged`, `kInvalid` is set.

**Not covered.**

- Relativistic mode beyond first order: `kPairwise` is non-relativistic, so a
  relativistic event gets `kBorn` (flag `kRelativistic`), with an error
  estimated from the non-relativistic second order; only `kBrute` (which
  integrates the relativistic copy) is exact there. Retardation is ignored,
  as in the event.
- Non-central forces (checked numerically at the present separation) and
  events with `SetPairSymmetric(kFALSE)`: first order only.
- Velocity-dependent forces, and forces that read particle data other than
  mass, charge and position (the pairs are evaluated with stand-in
  `EMparticle`s): not supported.
- Force types other than the four of `force.hpp` are handled by the generic
  quadrature and are untested; only `coulombForce` and `universalForce` of
  exactly that type (not subclasses) use the closed forms.
- `oddevent` is not modelled: the private copies are `event`s.
- Invincible particles are fixed (their stored velocity is ignored) and are
  left out of momentum and energy sums.

**How far must I go?** With `kPairwise` the error falls like \f$ 1/R^2 \f$
instead of \f$ 1/R \f$, and it is estimated. Stop with `ASYD` at
\f$ \epsilon \f$ (all pairs receding, remaining deflection below
\f$ \epsilon \f$): the pairwise error then scales like \f$ \epsilon^2 \f$.
Require `IsReliable()` and `GetMaxRelError()` below what you need, and check
a sample with `kBrute`.

# Forces and potentials {#con_forces}

Every force is central and pairwise, `Evaluate()` is exactly
\f$ -\nabla_{x_1} \f$ `Potential()` (checked in `test_force_check`), and
`Evaluate(p2, p1) = -Evaluate(p1, p2)`. Each force acts on every pair of the
event with the same parameters; one particle is one body.

| Class | Potential | Parameters (defaults) |
|---|---|---|
| \ref coulombForce | \f$ k q_1 q_2 / \sqrt{r^2+\varepsilon^2} \f$ | k = 0.00154629 amu c\f$^2\f$ fm/e\f$^2\f$ (`SetK`), softening \f$ \varepsilon \f$ = 0 fm (`SetSoftening`) |
| \ref universalForce | \f$ -G m_1 m_2 / \sqrt{r^2+\varepsilon^2} \f$ | G = 1.2331e-39 fm c\f$^2\f$/amu (`SetG`), \f$ \varepsilon \f$ = 0 |
| \ref yukawaForce | \f$ -g\, e^{-r/a}/r \f$ | g = 0.1 amu c\f$^2\f$ fm = 93.1 MeV fm (**illustrative, not fitted**), a = 1.4138 fm = \f$ \hbar/(m_\pi c) \f$ |
| \ref woodsSaxonForce | \f$ -V_0 / (1 + e^{(r-R)/a}) \f$ | V0 = 0.0536772 amu c\f$^2\f$ = 50 MeV, R = 8.886 fm = 1.2 (4^(1/3) + 197^(1/3)) fm, a = 0.65 fm |

- Coulomb reads the charges with `particle::GetQ()` (0 for a plain
  `particle`), so neutral particles feel no Coulomb force. With
  \f$ \varepsilon = 0 \f$ the original unsoftened expression is used, bit for
  bit; the softened force is finite at \f$ r = 0 \f$.
- Gravity uses rest masses and is always attractive. Newton's constant in
  these units: \f$ G m_u / c^2 = 1.2331\times10^{-54} \f$ m
  \f$ = 1.2331\times10^{-39} \f$ fm, so \f$ G \f$ = 1.2331e-39
  fm c\f$^2\f$/amu (the default; before November 2026 the code had
  1.2331e-17, a wrong exponent). At nuclear scales it is negligible.
- Yukawa: \f$ g > 0 \f$ attracts. It diverges at \f$ r \to 0 \f$ like Coulomb.
- Woods-Saxon: \f$ V_0 > 0 \f$ attracts; finite everywhere, the null vector at
  \f$ r = 0 \f$; written so that \f$ e^x \f$ never overflows.
- Your own force: see \ref tut_forces. It must be a function of positions and
  particle constants only.

# Path recording and closest approach {#con_path}

A particle created with `record = true` owns a \ref path, a list of
(position, velocity, time). After each step the state is appended if the
particle's step count is a multiple of its record period (`SetRecordPeriod(n)`,
default 1; values below 1 are taken as 1). The initial state is appended
only by an explicit `AddRecord()`, which appends regardless of the period.
`path::GetLast*()` is the last appended state; it equals the final state only
if the last step was recorded. An invincible recording particle records the
time 0 at every step.

`path::GetDCA(point, DCA, t)` finds the recorded sample with the smallest
squared distance \f$ d^2 \f$ to a fixed point and, by default, refines it
with the parabola through the samples \f$ k-1, k, k+1 \f$ (in \f$ d^2 \f$
against time, non-uniform times allowed): with
\f$ s_0 = t_{k-1} - t_k < 0 < s_2 = t_{k+1} - t_k \f$,
\f$ q_i = (d^2_{k\pm1} - d^2_k)/s_i \f$, \f$ a = (q_0 - q_2)/(s_0 - s_2) \f$,
\f$ b = q_2 - a s_2 \f$, the minimum is at \f$ t_k - b/(2a) \f$ with
\f$ d^2 = d^2_k - b^2/(4a) \f$. The discrete value is kept at the ends of the
path, for non-increasing times, for \f$ a \le 0 \f$ and for a vertex outside
\f$ [s_0, s_2] \f$. `GetDCA(path*, ...)` does the same for the distance
between two tracks, which must have the same length and identical times
(same record period, both moving). On failure both outputs are -1.

# Reproducibility {#con_repro}

A result file is reproducible when it records everything that determined it:

- **Seed.** `runinfo::SeedFromJob(jobIndex, seed = -1)` returns `seed` if it
  is given (\f$ \ge 0 \f$), otherwise `runinfo::kDefaultSeed` (65539) +
  jobIndex. Use one `TRandom3` per job, seeded with it.
- **Metadata.** `runinfo::WriteCommon(seed, dt, integrator, nevents,
  inspectors)` writes the `TNamed` keys `seed`, `dt`, `integrator`,
  `nevents`, `inspectors` and `git_sha` into the current ROOT directory;
  `runinfo::Write(key, value)` adds your own (string, double with 17
  significant digits, or integer). Write them before closing the file.
- **Code revision.** `git_sha` is the short commit hash CMake found when it
  configured the build (`RP_GIT_SHA`), `unknown` outside a git checkout.
  Uncommitted changes are not marked, and a new commit on the same branch is
  picked up only at the next CMake run.
- Integrator, `dt`, tolerances and inspectors change the numbers; write all of
  them (`WriteCommon` covers `dt`, the integrator and a text description of
  the inspectors; add tolerances and the stop distance yourself).

# Known limitations {#con_limits}

From `TODO.md` and from checking the code while writing this manual:

- `TAYLOR` is still the default integrator (decision pending).
- `macro/condor/condor_run.sh` passes `<title> $(ClusterId) $(ProcId)`, so
  the job index is the cluster id and the seed the process id: with
  `queue N` > 1 the processes of a cluster write the same file (\ref tut_mc).
- `macro/plot_rutherford.C` compares nominal `b` and `v`, so its residual
  shows the finite-distance bias, not the integration error.
- The Coulomb constant is \f$ 2.8\times10^{-4} \f$ high (section
  \ref con_units).
- An invincible particle's own time does not advance; `event::GetTime()` is
  the largest time over the particles.
- The Verlet force cache does not notice changes of force parameters
  (\ref con_integrators).
- `asymptote` (two-body) returns no asymptote for a receding radial pair
  under an attraction (`kRadialAttractive`); `asymptoteN` handles it.
- `inspectorP::SetInitial_D()` (prints that it needs fixing), `preEvent`
  (commented out) and `event::preDeriveSetup()` (empty) are unfinished; do
  not use them. \ref experiment and \ref experimentP (a template plus
  `MakeEvent(...)` with the new particles) work, but the test programs still
  use their own copies of such a helper.
- Inspector methods are strings parsed once into an enumerator; a typo is
  caught only when `DeriveInspect()` starts.
- The relativistic mode has no retardation or magnetic forces
  (\ref con_rel).
