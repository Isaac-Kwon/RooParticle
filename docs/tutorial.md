\page tutorial Tutorial

RooParticle integrates the classical motion of several charged bodies that
interact pairwise, until conditions that you choose say that the collision is
over. It is written for **three or more bodies**, where no formula gives the
answer. This tutorial starts with two bodies anyway, an alpha particle on a
gold-like nucleus, because that case has an exact answer (Rutherford) and lets
you see how accurate each step is. Section \ref tut_many then carries the same
tools to a fixed target with two projectiles and to three free bodies, and
section \ref tut_nbody_asymptote shows how to get their asymptotic velocities,
with an error estimate, without integrating to enormous distances.

Every code block on this page is a `\snippet` of a program in `test/` named
`main_tutorial_<topic>_check.cpp`. Each one is a complete program without
arguments that runs in well under a few seconds, prints what is shown here,
and exits with a non-zero code if its own results are wrong. They run in
`ctest`, so the code on this page compiles and works with the library it
documents. The printed output next to the code was copied from a run; your
last digits may differ with another compiler.

The background for each step (formulas, the full tables of options, what is
and is not conserved) is in \ref concepts "Concepts". Physics checks against
known results are in \ref validation "Validation".

Units everywhere: length fm, mass amu, velocity c, charge e, time fm/c,
energy amu c\f$^2\f$ (1 amu c\f$^2\f$ = 931.494 MeV). See \ref con_units.

# 0. Build, test, and build this manual {#tut_build}

You need ROOT, CMake and a C++ compiler (Linux or WSL). Doxygen is optional
and only needed for this manual.

~~~{.sh}
source <ROOT directory>/bin/thisroot.sh
mkdir build && cd build
cmake <RooParticle directory>
make -j8
ctest --output-on-failure     # the self-checking tests, including this tutorial
make doc                      # this manual: docs/html/index.html
~~~

`make` builds the library `librooparticle.a`, the program `main`, and one
executable `test_<name>` for every `test/main_<name>.cpp`, all in the build
directory. Every `test/main_*_check.cpp` (and `test_integrator_compare`) is
registered in `ctest` automatically when CMake runs; a new tutorial or check
program only needs a re-run of `cmake`. Programs write their ROOT files to
`Data/` below the working directory (CMake creates `build/Data`).

Doxygen warnings go to `build/docs/doxygen-warnings.log`, not to the terminal.

# 1. The first event {#tut_first}

Program: `test/main_tutorial_first_check.cpp`.

Three kinds of objects make an event:

- a **particle** (\ref particle) is a point mass with a position and a
  velocity. \ref EMparticle is a particle with a charge. A particle can be
  **invincible**: then it acts on the others but never moves (a fixed target,
  an infinitely heavy centre). A particle can **record** its track (a
  \ref path of position, velocity and time after each step).
- a **force** (\ref force) is a pairwise interaction. \ref coulombForce is
  \f$ U = k q_1 q_2 / r \f$.
- the **event** (\ref event) holds the particles and forces of one collision
  and advances them in time, one **step** of length `dt` at a time.

The headers:

\snippet main_tutorial_first_check.cpp includes

The setup. The event **owns** the particles you give to `AddParticle()`:
create them with `new` and let the event delete them. It does **not** own the
force: the force must live at least as long as the event, and it can be shared
by many events (see \ref con_ownership).

\snippet main_tutorial_first_check.cpp setup

The target mass (197 amu) plays no role here because the target is
invincible. Choose the time integrator. `VERLET` (velocity Verlet) is a good
first choice: one force evaluation per step and no energy drift. The default
is still `TAYLOR`, the original scheme, which is much less accurate (section
\ref tut_accuracy and \ref con_integrators).

\snippet main_tutorial_first_check.cpp integrator

Now derive (advance) the event for a fixed number of steps. The alpha needs
about \f$ 2 \times 1000\,\mathrm{fm} / 0.05\,c = 40000 \f$ fm/c to come in
and go out again to about 1000 fm. `DeriveDTN(n, dt)` makes `n` steps of
`dt`.

\snippet main_tutorial_first_check.cpp derive

Reading the final state: `GetX()`, `GetV()` (current position and velocity),
`GetECNT()` (number of steps this particle made), `GetTime()` (its evolved
time; invincible particles count steps but their time stays 0).

\snippet main_tutorial_first_check.cpp result

Energy conservation is the first thing to check after any run:
`event::GetNetE()` is the kinetic energy of the moving particles plus the
potential energy of all pairs.

\snippet main_tutorial_first_check.cpp energy

Output:

~~~{.txt}
integrator           : VERLET
steps, time          : 40000 steps, 40000.0 fm/c
final position (fm)  : (547.085, 753.648, 0.000)  r = 931.283 fm
final velocity (c)   : (0.031449, 0.038753, 0.000000)  |v| = 0.049908 c
scattering angle     : 50.9402 deg
energy  E0, E1       : 5.244008999145e-03  5.244008999121e-03 amu c^2
relative drift       : 4.54e-12
kinetic energy       : 4.6575 MeV at the start, 4.6404 MeV at the end
~~~

Two things to notice. The energy is conserved to \f$ 5 \times 10^{-12} \f$.
And the alpha ends at 931 fm, not 1000 fm, with a kinetic energy different
from the start: at 1000 fm the Coulomb energy is still about 5 % of the
kinetic energy, so neither the start nor the end is "far away". Section
\ref tut_physics deals with that. Choosing a fixed number of steps by hand
does not carry over to other impact parameters or to more bodies; the next
section lets the event decide when to stop.

# 2. Stopping automatically {#tut_stop}

Program: `test/main_tutorial_stop_check.cpp`.

An **inspector** is one terminating condition. `event::DeriveInspect()`
derives step after step and, before each step, asks the inspectors whether to
stop. There are two kinds:

- \ref inspectorP looks at a **pair** of particles: their relative position
  \f$ \Delta x = x_1 - x_2 \f$ and velocity \f$ \Delta v = v_1 - v_2 \f$
  (angle between them, distance, radial velocity, speed, step count).
- \ref inspectorE looks at the whole **event**: energies, time, smallest pair
  distance, largest acceleration, energy drift, leaving a volume.

Each inspector is a method name (a short string such as `"DEG"`) and a
criterion value. The full lists are in \ref con_stop. The event owns the
inspectors you add, like the particles.

All the events of this section use the setup of section 1 (alpha at
\f$ x = -1000 \f$ fm, impact parameter `b`, VERLET, `dt` = 1 fm/c) and print
one line with this function, where \f$ v_\infty = \sqrt{2E/m} \f$ is the speed
the alpha will have infinitely far away:

\snippet main_tutorial_stop_check.cpp report

## The default: AND {#tut_stop_and}

By default the inspectors are combined with **AND**: the event stops at the
first inspection where all of them are satisfied. This is the setup the
project has always used:

\snippet main_tutorial_stop_check.cpp and

`"DEG", 1` is satisfied when the angle between \f$ \Delta x \f$ and
\f$ \Delta v \f$ is below 1 degree, i.e. the alpha moves almost radially away
from the target: it has reached its asymptotic direction to about a degree.
`"CNT", 50000` is satisfied when both particles have made more than 50000
steps. `CNT` is a **minimum**, not a maximum.

Why AND, and why the minimum step count? Near a head-on collision the alpha
turns around and immediately moves straight out, so the angle condition is
satisfied right after the turning point, close to the target, where the
alpha is still slow:

\snippet main_tutorial_stop_check.cpp why

~~~{.txt}
b=50: DEG<1 AND CNT>50000          INSPECTOR    77436 steps  r =  2820.9 fm  |v|/v_inf = 0.99171
b=1:  DEG<1 alone                  INSPECTOR    22907 steps  r =    83.9 fm  |v|/v_inf = 0.66687
b=1:  DEG<1 AND CNT>50000          INSPECTOR    50001 steps  r =  1393.3 fm  |v|/v_inf = 0.98314
~~~

With `DEG<1` alone the b = 1 fm event stops 84 fm from the target with two
thirds of its final speed. The AND with `CNT>50000` makes sure every particle
travels far enough, whatever the impact parameter. That is the intent of the
default. Note also that even the b = 50 fm event, stopped at 2821 fm, is still
0.8 % below its final speed: the Coulomb tail is long (section
\ref tut_physics).

## An event inspector {#tut_stop_event}

The condition "the interaction energy has become small" belongs to the whole
event, so it is an \ref inspectorE. Here: the pair is receding (`"VRG", 0`:
radial velocity \f$ \Delta x \cdot \Delta v / |\Delta x| > 0 \f$) AND
\f$ |PE/KE| < 0.01 \f$:

\snippet main_tutorial_stop_check.cpp event

## OR mode {#tut_stop_or}

`SetInspectMode(event::kInspectOR)` stops at the first inspection where any
inspector is satisfied:

\snippet main_tutorial_stop_check.cpp or

## Abort inspectors {#tut_stop_abort}

An inspector marked with `SetAbort()` is an **abort inspector**: it stops the
event on its own, whatever the mode and the other inspectors say, and
`DeriveInspect()` then returns `kFALSE`. Use it for situations where the
result is not wanted or not trustworthy: a close encounter (`"DSL"`), a
numerical failure (`"EDG"`: relative energy drift above a value).

\snippet main_tutorial_stop_check.cpp abort

## The step cap, typos, and the stop reason {#tut_stop_cap}

`DeriveInspect()` never runs more than `SetMaxInspectSteps(n)` steps in one
call (default 100,000,000; 0 for no limit). Inspectors whose method name is
not defined are refused before anything is derived. A fifth way to end is
`kStopFailed`: the adaptive integrator `RK45` gave up inside a step (section
\ref tut_forces, softening).

\snippet main_tutorial_stop_check.cpp cap

\snippet main_tutorial_stop_check.cpp invalid

`event::GetStopReason()` says why the last `DeriveInspect()` returned. In a
production loop, act on it:

\snippet main_tutorial_stop_check.cpp switch

Output of the remaining events (the cap and the typo also print a message
on stderr):

~~~{.txt}
b=50: VRG>0 AND PKL<0.01           INSPECTOR   114480 steps  r =  4705.5 fm  |v|/v_inf = 0.99504
b=50: DSTG>3000 OR TMG>30000       INSPECTOR    30001 steps  r =   439.4 fm  |v|/v_inf = 0.94551
b=1:  DSTG>3000, abort DSL<60      ABORT        20063 steps  r =    60.0 fm  |v|/v_inf = 0.47251
b=50: DSTG>3000, cap 10000         MAXSTEPS     10000 steps  r =   507.2 fm  |v|/v_inf = 0.95296
b=50: "DGE" (typo)                 INVALID          0 steps  r =  1001.2 fm  |v|/v_inf = 0.97646
b=50: same, inspected every 10     INSPECTOR    77440 steps  r =  2821.1 fm  |v|/v_inf = 0.99171
~~~

Things that surprise people:

- The arguments are `DeriveInspect(iperiod, dt)`: the **inspection period
  first**, the time step second, both with default 1. `DeriveInspect(0.5)`
  does not use a step of 0.5 fm/c.
- The inspectors are evaluated before the first step too. If they are
  already satisfied, nothing is derived. An event with no inspectors at all
  stops at once; an event with only abort inspectors runs until one fires
  (or the cap).
- `CNT` counts the steps each particle has made since it was created, and a
  copied particle keeps its count (section \ref tut_mc). The cap counts the
  steps of the current call only.
- The angle methods `DEG`, `RAD`, `COS` are never satisfied while
  \f$ \Delta x \f$ or \f$ \Delta v \f$ is exactly zero (two particles at the
  same place, or with the same velocity, e.g. both at rest), and `VRG`/`VRL`
  never while \f$ \Delta x = 0 \f$.
- `"VOUT"`/`"VOUTA"` (leaving the event volume) are never satisfied when no
  volume was given with `SetVolume()`.

# 3. Getting the physics out (two bodies) {#tut_physics}

Program: `test/main_tutorial_physics_check.cpp`. Same alpha, b = 50 fm,
stopped by `DEG<1 AND CNT>50000`, now with `record = true` so that its track
is kept.

\snippet main_tutorial_physics_check.cpp run

## The raw scattering angle

The angle between the initial and the final velocity:

\snippet main_tutorial_physics_check.cpp raw

~~~{.txt}
raw angle at r = 2821 fm          : 50.975599 deg
~~~

## Why the nominal b and v are not the asymptotic ones {#tut_physics_finite}

The Rutherford formula
\f[ \tan\frac{\theta}{2} = \frac{\kappa}{\mu\, v_\infty^2\, b}, \qquad \kappa = k q_1 q_2 \f]
holds for the speed \f$ v_\infty \f$ and impact parameter \f$ b \f$ of the
**asymptotes**, i.e. infinitely far away. The simulation starts at 1000 fm
with \f$ v = 0.05\,c \f$ and \f$ y = 50 \f$ fm. There the potential energy is
about 5 % of the kinetic energy, so the speed at infinity is larger, and the
incoming straight line the particle came on (its incoming asymptote) is not
the line \f$ y = 50 \f$ fm. Likewise the final velocity at 2821 fm is not yet
the asymptotic one. The nominal values are only approximations.

## The exact two-body answer: asymptote {#tut_physics_asymptote}

\ref asymptote takes the state of **two** bodies at **any one instant** and
returns the exact incoming and outgoing asymptotes, using the conserved
energy, angular momentum and Laplace-Runge-Lenz vector of the pair (formulas
in \ref con_asymptote). `kappa` is the strength of \f$ U = \kappa / r \f$:

\snippet main_tutorial_physics_check.cpp asymptote

~~~{.txt}
status                          : UNBOUND
thetaCM, extrapolated at stop   : 51.013865 deg
thetaCM, from the initial state : 51.013865 deg
v_inf                           : 0.0512055 c   (nominal 0.05 c)
effective b = |L| / (mu v_inf)  : 48.8229 fm  (nominal 50.0 fm)
tilt of the incoming asymptote  : 0.0341 deg
lab angle vs the x axis at +inf : 50.979763 deg
Rutherford at nominal (b, v)    : 52.082826 deg
Rutherford at (b_eff, v_inf)    : 51.013865 deg
~~~

Read this carefully:

- The asymptote computed before any step and the one computed at the stop
  agree: for two bodies you do not need to integrate at all. Integrating and
  comparing is a test of the integrator.
- The nominal values give the Rutherford angle 1.07 degree too large. With
  the effective \f$ b \f$ and \f$ v_\infty \f$ the formula is exact.
- `GetThetaCM()` is the angle between the incoming and the outgoing
  asymptote. The incoming asymptote is tilted by 0.034 degree from the x
  axis, so the angle measured from the initial velocity at +infinity
  (50.980 degree) differs from it by that tilt; the raw angle at 2821 fm
  (50.976 degree) has not quite reached that value.

\warning `asymptote` is **exact only for exactly two bodies with a pure
\f$ 1/r \f$ interaction**, non-relativistic. A third body, a nuclear or
screened force, softening, or the relativistic mode all break it. For a
short-range extra force it is still exact once that force has become
negligible (evaluate it far from the nucleus, not at the start). With three
or more bodies do not use it on a pair: use \ref asymptoteN, which takes the
whole event and returns every asymptotic velocity with an error estimate
(section \ref tut_nbody_asymptote). `asymptote` also returns no asymptote for
a receding pair on a straight radial line under an attraction (status
`kRadialAttractive`), although one exists; `asymptoteN` handles that case.

## Distance of closest approach

The recorded \ref path gives the distance of closest approach to a point
(here the fixed target) with `path::GetDCA()`: the smallest distance over the
recorded steps, refined by a parabola through the three samples around it.

\snippet main_tutorial_physics_check.cpp dca

~~~{.txt}
DCA, discrete                   : 77.389937188 fm at t = 20701.0 fm/c
DCA, interpolated               : 77.389936549 fm at t = 20701.232 fm/c
closest approach, asymptote     : 77.389937216 fm
path points                     : 77436
~~~

The interpolated value is the minimum of the computed trajectory between the
samples; it differs from the exact closest approach by the integration error
of VERLET at `dt` = 1 (\f$ 7 \times 10^{-7} \f$ fm). With RK4 the two agree
to \f$ 10^{-12} \f$ fm (section \ref tut_record). For two moving particles
there is `GetDCA(path*)`, which needs both paths recorded on the same steps.

## Starting at the nominal energy: CalibrateVelocity()

`event::CalibrateVelocity()` rescales the speed of every moving particle
(keeping its direction) so that the total energy equals the kinetic energy you
gave. Then the nominal speed is the speed at infinity:

\snippet main_tutorial_physics_check.cpp calibrate

~~~{.txt}
after CalibrateVelocity         : |v| = 0.0487647 c, v_inf = 0.050000000 c, b_eff = 48.7647 fm
~~~

The impact parameter is still not the nominal one. When several particles
move, the potential energy of a pair is shared between them, which keeps the
total energy right but not each particle's own asymptotic speed. A start
inside the Coulomb barrier (potential energy larger than the kinetic energy)
gives NaN.

# 4. Three and more bodies {#tut_many}

Program: `test/main_tutorial_multibody_check.cpp`. This is what the code is
for. Three things change when more than two bodies interact:

1. There is no closed formula. Conservation laws are the correctness
   evidence inside one run; convergence (in `dt`, tolerance, stopping
   distance) or the error estimate of \ref asymptoteN is the evidence for the
   numbers you keep.
2. Terminating conditions must cover every pair, or the whole event.
3. The update order matters: `oddevent` gives different results from `event`.

The program of this section shows the events and the brute way to their
asymptotic velocities (integrate far). The next section,
\ref tut_nbody_asymptote, shows the way to use in practice.

## Terminating many bodies

A pair inspector only knows about its two particles. `"DEG"` on one pair
says nothing about the others, and with several projectiles the "angle
between position and velocity" of a pair is not an asymptote criterion at
all. What works for any number of bodies:

- the event inspector `"ASYD"` (or `"ASYV"`): every moving particle will
  still turn by less than the value (rad), or change its velocity by less
  than that fraction, from now to infinity, **and** every interacting pair is
  receding. This is the condition meant for three and more bodies; it is
  used in section \ref tut_nbody_asymptote.
- one `"VRG", 0` (receding) inspector for every pair, combined with AND,
- event inspectors: `"PKL"` (interaction energy small compared to kinetic
  energy), `"ACL"` (largest acceleration small), `"DSG"` (smallest pair
  distance large), `"VOUTA"` (all moving particles out of a box),
- an abort on `"EDG"` (energy drift: a numerical failure) and a step cap.

This section uses the second group, which says nothing about how close the
velocities are to their final values:

\snippet main_tutorial_multibody_check.cpp pairs

\snippet main_tutorial_multibody_check.cpp termination

`"DSG"` and `"PKL"` include every pair, also pairs of two invincible
particles: with two fixed centres, `DSG` cannot exceed their separation.

## A. Fixed target and two projectiles

Two alpha particles come in together, 70 fm apart, and repel each other as
well as the fixed gold-like target.

\snippet main_tutorial_multibody_check.cpp fixed

~~~{.txt}
inspectors          : 5 (3 VRG pairs, PKL, abort EDG)
stop                : INSPECTOR after 115268 steps
|dE/E|              : 1.62e-11
|dL|/|L| (about Au) : 1.44e-12
|P| (not conserved, the fixed centre takes the recoil) : 0.4000 -> 0.2833 amu c
~~~

With an invincible particle the total momentum is **not** conserved (the
fixed centre absorbs momentum without moving), but the angular momentum
**about the fixed centre** is, because every force is central: pass the
centre's position to `GetNetAngularMomentum()`. Energy is always conserved.

The effect of the alpha-alpha force, from the same alphas run alone with the
target:

\snippet main_tutorial_multibody_check.cpp alone

~~~{.txt}
alpha 1 : angle 41.5327 deg with the other alpha, 61.5951 deg alone, |v| 0.052842 vs 0.050952 c
alpha 2 : angle 51.9637 deg with the other alpha, 78.6386 deg alone, |v| 0.049267 vs 0.050885 c
~~~

At this distance the alphas change each other's angles by 20 to 27 degrees
and exchange energy; a two-body treatment of each alpha would be wrong.

## B. All bodies free, and asymptotic velocities by integrating far

A 12C nucleus at rest, an alpha and a proton coming in; nothing is fixed.
Use the adaptive integrator `RK45` (section \ref tut_accuracy): it takes small
substeps during the encounter and large ones far away. Energy, momentum and
angular momentum are all conserved now.

\snippet main_tutorial_multibody_check.cpp free

To get the asymptotic velocities, integrate until all pairs recede and the
closest pair is `R` apart, record the velocities, raise `R` tenfold and
continue the same event. The change from one `R` to the next falls like
\f$ 1/R \f$ (the Coulomb tail), which is the convergence check;
\f$ |PE/E| \f$ tells how much energy is still in the interaction.

\snippet main_tutorial_multibody_check.cpp converge

~~~{.txt}
   R (fm)    |v C| (c)    |v alpha|        |v p|    angle a    angle p     |PE/E|
    1e+03  0.001604310  0.049856513  0.041219732   8.202256  13.176939   2.83e-03   max rel. change 0.0e+00
    1e+04  0.001601719  0.049916597  0.041296944   8.237264  13.371705   2.51e-04   max rel. change 1.4e-02
    1e+05  0.001601535  0.049922033  0.041302762   8.240715  13.390723   2.48e-05   max rel. change 1.2e-03
    1e+06  0.001601517  0.049922572  0.041303331   8.241059  13.392622   2.48e-06   max rel. change 1.1e-04
    1e+07  0.001601515  0.049922626  0.041303387   8.241094  13.392812   2.48e-07   max rel. change 1.1e-05
extrapolated |v_inf| of particle 0 : 0.001601515 c (correction 1.3e-06)
extrapolated |v_inf| of particle 1 : 0.049922632 c (correction 1.4e-07)
extrapolated |v_inf| of particle 2 : 0.041303394 c (correction 4.0e-07)
~~~

Stopping at 1000 fm would give the proton direction wrong by 0.2 degree;
\f$ 10^7 \f$ fm brings the change below \f$ 2 \times 10^{-5} \f$, and the
\f$ 1/R \f$ (Richardson) extrapolation removes most of what is left. The
conservation laws at the end:

\snippet main_tutorial_multibody_check.cpp conserved

~~~{.txt}
|dE/E|              : 8.15e-15
|dP| / sum |p|      : 2.44e-15
|dL| / |L|          : 1.24e-09
RK45                : 25417 steps, 25421 substeps, 0 rejected
~~~

They say the integration did what it was asked; they do not say the stopping
distance was large enough. That is what the table above is for. When you
compare two runs (two tolerances, two integrators), compare them at the same
stopping point: an overshoot of the stop by 100 fm at 1e4 fm already moves
the alpha direction by about \f$ 10^{-5} \f$ degree.

Integrating to \f$ 10^7 \f$ fm works, but it is slow and the accuracy is only
known from the table. The next section gets the same answer at 1400 fm.

## Asymptotic velocities of many bodies: asymptoteN {#tut_nbody_asymptote}

Program: `test/main_tutorial_asymptoteN_check.cpp`, the three free bodies of
part B.

\ref asymptoteN takes the present state of a whole event (any number of
bodies, any of the forces of `force.hpp`) and returns, for every moving
particle, the asymptotic velocity \f$ v_\infty \f$ **and an estimate of its
error**. The event is not changed. The default method, `kPairwise`, solves
every pair exactly as an isolated two-body problem and adds up the velocity
changes; what it misses is the three-body cross term, which falls like
\f$ 1/R^2 \f$ with the distance \f$ R \f$ at which you stop
(\ref con_nbody_asymptote). It needs no integration: a few microseconds per
call. The method `kBrute` continues a private copy of the event with RK45 in
doubling time blocks and extrapolates; it is the independent reference used
below. Use `asymptoteN` for three or more bodies, and \ref asymptote only for
the exact two-body case.

### Stop with ASYD, read with asymptoteN

Stop when every particle has less than a given deflection left (inspector
`"ASYD"`), then read the asymptotes. The program does this for three values
of `ASYD` and compares each result with `kBrute` from the same state:

\snippet main_tutorial_asymptoteN_check.cpp asyd

~~~{.txt}
   ASYD stop (fm)    steps   sweeps    true turn  raw rel err pair rel err  estimate reliable
  1e-02      1414      873     6405     9.98e-03      1.0e-02      3.7e-06   2.4e-05      yes
  1e-03     12768     6832    48118     1.00e-03      1.0e-03      4.0e-08   2.6e-07      yes
  1e-04    126464    66401   465101     1.00e-04      1.0e-04      4.0e-10   2.6e-09      yes
~~~

Columns: closest pair distance at the stop; steps and force sweeps (one
sweep = all pairs evaluated once) of the run; the true remaining deflection
(angle between the velocity at the stop and the brute-force \f$ v_\infty \f$);
the largest relative error of the raw velocities and of the pairwise
\f$ v_\infty \f$ against the brute force; the pairwise error estimate
(`GetMaxRelError()`); `IsReliable()`.

- `ASYD` means what it says: the true remaining deflection is the criterion
  value.
- The raw velocities at the stop are off by about the same relative amount.
  The pairwise correction removes all but \f$ 4 \times 10^{-6} \f$ of it at
  1414 fm, and its error falls like \f$ 1/R^2 \f$ (100 times per factor 10 in
  distance) while the raw error falls like \f$ 1/R \f$.
- The estimate is larger than the actual error by a factor 6 to 7 here
  (`asymptoteN.hpp` reports 4 to 9 in its own tests): use it as an upper
  bound, not as the error itself.
- The relative criteria are per particle. The slow recoiling 12C
  (\f$ 0.0016\,c \f$) dominates `ASYD`: its direction still changes by
  \f$ 10^{-2} \f$ rad at 1414 fm, while the alpha's changes by
  \f$ 5 \times 10^{-4} \f$.

Details of the stop at `ASYD` = \f$ 10^{-2} \f$:

\snippet main_tutorial_asymptoteN_check.cpp pairwise

~~~{.txt}
stop                : INSPECTOR, closest pair 1414 fm, remaining deflection (first order) 9.992e-03 rad
pairwise            : reliable yes, flags [OK], largest relative error 2.4e-05
  particle 0 : |v| now 0.001603283 c, |v_inf| 0.001601519 c +- 5.0e-08, still to turn 9.99e-03 rad, pairwise - brute 4.5e-09 c
  particle 1 : |v| now 0.049877347 c, |v_inf| 0.049922655 c +- 2.2e-07, still to turn 4.77e-04 rad, pairwise - brute 3.1e-08 c
  particle 2 : |v| now 0.041249346 c, |v_inf| 0.041303279 c +- 1.2e-06, still to turn 2.64e-03 rad, pairwise - brute 1.5e-07 c
brute force         : 212 RK45 substeps (1286 sweeps), integrated 1.1e+09 fm/c on
~~~

`"ASYV"` is the same with the relative velocity change
\f$ \max_i |v_{\infty,i} - v_i| / |v_{\infty,i}| \f$:

\snippet main_tutorial_asymptoteN_check.cpp asyv

~~~{.txt}
ASYV < 1e-3         : INSPECTOR at 12830 fm, true remaining relative change 1.00e-03
~~~

### What it saves

From the `ASYD` = \f$ 10^{-2} \f$ stop, the program continues a copy of the
event the way part B does (calls of `DeriveDT(R/100)` up to distance `R`) and
compares the raw velocities and their \f$ 1/R \f$ extrapolation with the
brute-force answer:

\snippet main_tutorial_asymptoteN_check.cpp far

~~~{.txt}
pairwise at the stop: relative error 3.7e-06, 0 force sweeps (closed form)
   R (fm)     sweeps    raw rel err     Richardson
    1e+04      31563        1.3e-03       -1.0e+00
    1e+05      64581        1.3e-04        1.8e-06
    1e+06      97592        1.3e-05        1.5e-08
    1e+07     130603        1.3e-06        1.6e-10
    1e+08     163614        1.3e-07        1.5e-12
    1e+09     196625        1.3e-08        1.4e-12
pairwise accuracy reached: extrapolated at R = 1e+05 fm after 64581 more sweeps, raw at R = 1e+07 fm after 130603 more sweeps
~~~

For the accuracy that `kPairwise` gives at 1414 fm (\f$ 4 \times 10^{-6} \f$),
the run must go on to \f$ 10^5 \f$ fm with the extrapolation (64581 more force
sweeps, 10 times the 6405 of the whole run up to the stop) or to
\f$ 10^7 \f$ fm without it (130603 more, 20 times). The `kBrute` reference
gets \f$ 10^{-11} \f$ in 1286 sweeps because it doubles its time blocks;
that is the efficient form of "integrate far" when you need a reference.

Guidance: stop with `ASYD` (or `ASYV`) at a value whose square is below the
accuracy you need, read `asymptoteN` (pairwise), and require `IsReliable()`
and `GetMaxRelError()` below your tolerance. Check a sample of events with
`kBrute`.

### Incoming asymptotes

A run started at 1000 fm does not start on its incoming asymptotes either.
`kIncoming` gives the true velocities at \f$ t \to -\infty \f$ (by time
reversal), from the initial state:

\snippet main_tutorial_asymptoteN_check.cpp incoming

~~~{.txt}
  particle 1 : nominal |v| 0.050000 c, incoming |v| 0.050615744 c, tilt 0.2378 deg, pairwise - brute 2.2e-06 c (estimate 1.4e-05)
  particle 2 : nominal |v| 0.040000 c, incoming |v| 0.038174940 c, tilt 1.2687 deg, pairwise - brute 6.4e-06 c (estimate 7.7e-05)
  12C        : at rest at the start, incoming |v| 5.22e-05 c
~~~

The proton, which starts 100 fm ahead of the alpha, was 4.6 % slower and came in 1.3
degree off the x axis; the "12C at rest" was moving. Quote these, not the
nominal values, as the initial conditions of the collision. At the start the
alpha and the proton are only 128 fm apart, so the incoming result is less
accurate than the outgoing one: actual error \f$ 2 \times 10^{-4} \f$
relative, estimate \f$ 2 \times 10^{-3} \f$. Start farther apart if you need
better.

### Flags

`GetFlags()` (and `GetFlagNames()`) say how a result was obtained and when
not to trust it. `IsReliable()` is false if any of `STRONG`, `NOT_ISOLATED`,
`MARGINAL`, `SINGULAR`, `NOT_CONVERGED` or `INVALID` is set.

\snippet main_tutorial_asymptoteN_check.cpp flags

~~~{.txt}
outgoing from t = 0 : separating no, flags [APPROACHING|STRONG], reliable no, largest error 7.5e-03 c
bound pair          : flags [APPROACHING|CLUSTER] [APPROACHING|CLUSTER] [APPROACHING], clusters 1, cluster |V_inf| 0.020000 c, reliable yes
bound to fixed      : light particle [CLUSTER|BOUND_TO_FIXED], v_inf 0.0 c ; alpha [APPROACHING]
~~~

- Before the collision (`APPROACHING`), the pairwise solution still predicts
  an outgoing state, but the other bodies will disturb the pairs strongly
  (`STRONG`): not reliable. Integrate through the collision first; `ASYD`
  waits until all pairs recede.
- A bound pair is merged into a cluster (`CLUSTER`): both members get the
  cluster's centre-of-mass velocity, and its internal motion is not an
  asymptote. `NOT_ISOLATED` would mean the other bodies may still break it.
- A body bound to an invincible particle has no asymptote (`BOUND_TO_FIXED`,
  \f$ v_\infty = 0 \f$).
- Other flags: `RELATIVISTIC` (first order only), `NO_PAIRWISE` (a non-central
  force or `SetPairSymmetric(kFALSE)`: first order only), `MARGINAL` (a pair
  with energy near zero), `SINGULAR` (radial fall into an attractive centre).

Two practical points. `asymptoteN` recognises the plain \ref coulombForce
and \ref universalForce by their exact type; any other force, **including
your own subclass of `coulombForce`**, is handled by numerical quadrature,
about 60 times slower per call (the `ASYD` inspection is evaluated at every
inspection). And the relativistic mode gets only the first-order (`kBorn`)
correction.

## C. Why oddevent differs once more than one body moves

\ref oddevent is a variant of the event that moves each particle **right
after** its own force is computed, so the later particles see the earlier
ones already moved. It has its own first-order scheme and ignores
`SetIntegrator()`, `SetPairSymmetric()` and `SetRelativistic()`.

\snippet main_tutorial_multibody_check.cpp odd

~~~{.txt}
one moving body     : |x_event - x_oddevent| = 0.0e+00 fm
three free bodies   : |dP| event 1.8e-15, oddevent 2.3e-06 amu c ; |x_event - x_oddevent| (alpha) = 8.10e-03 fm
~~~

With one moving body on a fixed target there is nobody to see a moved
particle, and `oddevent` is bit for bit `event` with `TAYLOR`. With several
moving bodies the action and the reaction of a pair are evaluated at
different states, so the total momentum drifts (here \f$ 2 \times 10^{-6} \f$
amu c over 60000 steps) and the trajectories differ. For many-body work use
`event` with a modern integrator; `oddevent` is kept because the program
`main` and older results use it.

# 5. Many events {#tut_mc}

Program: `test/main_tutorial_montecarlo_check.cpp`. A Monte Carlo over the
impact parameter for alpha on gold, written to a `TTree`. The same structure
serves any multi-body study: only the template and the per-event particles
change.

## Template event and copies

A **template event** holds what all events share. It is never derived; each
event is a copy (`new event(templ)`), to which the per-event particles and
inspectors are added.

\snippet main_tutorial_montecarlo_check.cpp template

Copying an event copies its particles (as new objects, with their state,
step count, charge and track), its inspectors (re-linked to the copied
particles) and its settings (integrator, tolerances, inspect mode, cap, pair
symmetry, relativistic flag). It does **not** copy the forces and the volume:
the copy uses the same objects. Two traps (details in \ref con_copy):

- Every particle is copied as a plain \ref particle. The charge survives
  (it lives in `particle`), but your own particle subclass would lose its
  type and data.
- `new event(*oddTemplate)` makes an `event`, not an `oddevent`. Copy an
  `oddevent` with `new oddevent(*oddTemplate)`.

## Seeds and job index

\snippet main_tutorial_montecarlo_check.cpp args

`runinfo::SeedFromJob(jobIndex)` gives 65539 + jobIndex unless you pass an
explicit seed, so parallel jobs never share a seed and rerunning a job index
reproduces it.

## The event loop and the TTree

\snippet main_tutorial_montecarlo_check.cpp tree

The impact parameter is drawn uniformly in the beam area
(\f$ b = b_{max} \sqrt{u} \f$), which is what a beam does; drawing `b`
uniformly (as `main` does) oversamples small `b` and must be weighted by
`b` for a cross section. Each event stops when the alpha recedes and is back
at 1000 fm; the exact angle then comes from `asymptote`, and the raw one is
stored next to it.

\snippet main_tutorial_montecarlo_check.cpp loop

`delete ev` deletes the event's particles and inspectors (and nothing else).

## Run metadata

`runinfo` writes the run parameters into the ROOT file as `TNamed` objects
(name = key, title = value), next to the tree, so the file says how it was
made. `WriteCommon()` writes `seed`, `dt`, `integrator`, `nevents`,
`inspectors` and `git_sha` (the revision the binary was built from;
`unknown` when the source is not a git checkout).

\snippet main_tutorial_montecarlo_check.cpp meta

\snippet main_tutorial_montecarlo_check.cpp readback

~~~{.txt}
file Data/tutorial_0.root : 200 events, seed 65539, git unknown
************************************************************************************************
*    Row   *         imp *        bEff *      SAngle *    thetaAsy *        rMin *       nstep *
************************************************************************************************
*        0 *   170.01774 *   166.06499 *  0.27486459 *  0.27889436 *   191.00136 *       40114 *
*        1 *   112.80338 *   110.16067 *  0.41419253 *  0.41687542 *   135.89785 *       40621 *
*        2 *   173.90557 *   169.86496 *  0.26862036 *   0.2727413 *   194.76604 *       40081 *
*        3 *   109.50609 *   106.93972 *  0.42644312 *  0.42904805 *   132.74818 *       40653 *
*        4 *   72.472261 *   70.768329 *   0.6343028 *  0.63602937 *    97.80016 *       41062 *
************************************************************************************************
events not ended by the inspectors : 0
max |thetaAsy - Rutherford(bEff, vInf)| : 2.2e-16 rad
~~~

200 events take about 2 s.

## Running jobs in parallel: macro/runmulti.py

`macro/runmulti.py target title njob [ncore]` runs
`target title <job index>` for job index 0 to `njob`-1, at most `ncore` at a
time (default 4), and exits non-zero if any job failed. Extra arguments go
after `--` and are put between the title and the job index.

~~~{.sh}
cd build
python3 <RooParticle directory>/macro/runmulti.py ./test_tutorial_montecarlo_check mc 4 -j 4
~~~

~~~{.txt}
RUN  job 0 (1 running)
RUN  job 1 (2 running)
RUN  job 2 (3 running)
RUN  job 3 (4 running)
DONE job 0 (exit 0, 3 running)
DONE job 1 (exit 0, 2 running)
DONE job 2 (exit 0, 1 running)
DONE job 3 (exit 0, 0 running)
all 4 jobs finished
~~~

This writes `Data/mc_0.root` ... `Data/mc_3.root` with seeds 65539 ... 65542.
Your program must take the job index as its **last** argument, as this one
does.

The programs `main` and `test_collider_HeHe` take
`[title] [jobIndex [seed]]`: a first argument that is not a number is the
title (it names the output file `Data/<title>_<jobIndex>.root`), so both
`./main 3` and `runmulti.py ./main <title> N` work. The Condor wrapper
`macro/condor/condor_run.sh` passes `<title> $(ClusterId) $(ProcId)`, i.e.
job index = cluster id and seed = process id: with `queue N` larger than 1
all processes of a cluster would write the same file.

To merge the jobs, `hadd Data/mc_all.root Data/mc_[0-3].root` works for the
trees. The metadata are kept as cycles (`seed;1` ... `seed;4`), and
`Get("seed")` returns only the last one; keep the job files if you need
each job's parameters.

## Plotting: macro/plot_rutherford.C

`macro/plot_rutherford.C` draws the branch `SAngle` [rad] against `imp` [fm]
with the Rutherford curve at the nominal speed, and the residual:

~~~{.sh}
root -l -b -q '<RooParticle directory>/macro/plot_rutherford.C("Data/mc_all.root", "rutherford", "Data/mc_all.png")'
~~~

~~~{.txt}
plot_rutherford: 800 events, v=0.05 m=4 q1=2 q2=79
  max |residual| = 0.024422 rad (signed -0.024422), RMS residual = 0.013045 rad
  saved Data/mc_all.png
~~~

The residual of up to 0.024 rad (1.4 degree) is not integration error: it is
the difference between nominal and asymptotic `b` and `v` explained in
section \ref tut_physics. The macro reads only the branches `imp` and
`SAngle`; to compare the exact values, draw `thetaAsy` against `bEff`
yourself (`tree->Draw("thetaAsy:bEff")`).

# 6. Choosing the accuracy {#tut_accuracy}

Program: `test/main_tutorial_accuracy_check.cpp`.

## Fixed step against adaptive

Every integrator except `RK45` uses the step `dt` you give. `RK45`
(Dormand-Prince 5(4)) still advances exactly `dt` and counts one step per
`DeriveDT(dt)`, but inside it takes as many substeps as its tolerance needs,
set with `SetTolerance(rel, abs)` (defaults \f$ 10^{-8} \f$ and
\f$ 10^{-10} \f$, see \ref con_rk45).

The cost of a run is the number of force evaluations. This force subclass
counts them:

\snippet main_tutorial_accuracy_check.cpp counting

Alpha on fixed gold, b = 50 fm, until receding at 1000 fm. The error is the
change of the asymptotic angle between the initial and the final state, i.e.
how well the run preserves the orbit:

\snippet main_tutorial_accuracy_check.cpp table

~~~{.txt}
integ       dt  rel tol     steps force eval     |dE/E| angle err (deg)
TAYLOR     1.0    0e+00     41377      41377    1.1e-04         7.9e-03
TAYLOR     0.1    0e+00    413774     413774    1.1e-05         7.9e-04
SEULER     1.0    0e+00     41378      41378    2.3e-06         5.2e-05
SEULER     0.1    0e+00    413775     413775    2.3e-07         5.2e-06
VERLET     4.0    0e+00     10345      10346    1.0e-12         2.1e-11
VERLET     1.0    0e+00     41378      41379    7.0e-14         1.5e-12
RK4        8.0    0e+00      5173      20692    3.3e-16         7.1e-13
RK4        2.0    0e+00     20689      82756    5.6e-15         1.8e-12
RK45     100.0    1e-06       414       2898    5.9e-11         1.9e-09
RK45     100.0    1e-10       414       3090    1.7e-11         4.9e-10
RK45    2000.0    1e-06        21        243    6.5e-07         1.8e-05
RK45    2000.0    1e-08        21        453    6.9e-09         1.7e-07
RK45    2000.0    1e-10        21        879    1.3e-10         3.1e-09
~~~

What to read from it:

- `TAYLOR` (the default) and `SEULER` are first order: ten times smaller
  `dt`, ten times smaller error, ten times the cost.
- `VERLET` costs one evaluation per step like them and is second order and
  symplectic (no energy drift). It is the sensible choice for fixed steps.
- `RK4` costs four evaluations per step and is fourth order.
- `RK45` with `dt` = 100 fm/c makes one substep per call: the **call** `dt`
  limits the substep, and the tolerance hardly matters. Give it a `dt` larger
  than the substep it would choose (here 2000 fm/c); then the tolerance
  controls the error, and the whole event costs a few hundred evaluations.
  Under `DeriveInspect()` the `dt` is also the inspection granularity: the
  stop can overshoot by up to one `dt` (here 100 fm at 0.05 c).

The measure matters. `test_integrator_compare` measures the error of the
direction at a fixed time instead (40000 fm/c), which includes the error of
the phase along the orbit; there `VERLET` is \f$ 3 \times 10^{-5} \f$ degree
at `dt` = 8 while `RK4` is \f$ 7 \times 10^{-11} \f$ degree (table in
\ref con_integrators). VERLET keeps the orbit (energy, angular momentum) far
better than its position along the orbit.

## No exact answer: convergence

For three free bodies (part B of section \ref tut_many) the same run at
tighter and tighter tolerance shows the error by its convergence:

\snippet main_tutorial_accuracy_check.cpp nbody

\snippet main_tutorial_accuracy_check.cpp nbodytable

~~~{.txt}
 rel tol   force eval alpha angle (deg)       change (deg)
   1e-05         5904     8.2371481169            0.0e+00
   1e-07         6552     8.2372775827            1.3e-04
   1e-09         7596     8.2372773564            2.3e-07
   1e-11        10872     8.2372773526            3.8e-09
   1e-13        19584     8.2372773526            4.4e-11
~~~

The angle at 1e4 fm is known to about \f$ 10^{-10} \f$ degree at
`rel` = \f$ 10^{-11} \f$; the stopping distance, not the integrator, is then
the limit (section \ref tut_many, part B).

## When oddevent differs

`oddevent` gives the same numbers as `event` with `TAYLOR` when only one
particle moves, and different, momentum-violating numbers when several move
(section \ref tut_many, part C). It is first order like `TAYLOR`.

# 7. Beyond a fixed Coulomb centre {#tut_forces}

Program: `test/main_tutorial_forces_check.cpp`.

## A recoiling target

Make the target an ordinary (not invincible) particle. `asymptote` then works
in the centre-of-mass frame and also gives the lab angles of both bodies:

\snippet main_tutorial_forces_check.cpp recoil

~~~{.txt}
fixed 12C : thetaCM 10.5817 deg, alpha lab 10.5817 deg, 12C lab 0.0000 deg, |dP|/|P| 1.8e-01
free 12C  : thetaCM 14.0693 deg, alpha lab 10.5718 deg, 12C lab 97.0344 deg, |dP|/|P| 3.5e-15
~~~

With a free 12C the momentum is conserved and the centre-of-mass angle is
larger (reduced mass 3 amu instead of 4). The 12C lab angle needs care. The
12C is at rest at the start, but the alpha already repels it there, so on its
incoming asymptote (\f$ t \to -\infty \f$) it was moving slowly toward the
alpha. `GetThetaLab2()` then measures from that incoming direction of its
own and gives 97.03 degree = 180 - 82.97 degree. The textbook recoil angle
\f$ (\pi - \theta_{CM})/2 = 82.97 \f$ degree (`asymptote::RecoilAngle()`)
assumes a target at rest at infinity. To simulate that, shift all velocities
by the target's asymptotic incoming velocity `GetV2In()` before deriving, as
`test/main_asymptote_check.cpp` does.

## A short-range nuclear force

Add a second force to the event; the forces are summed. At
\f$ v = 0.12\,c \f$ (26.8 MeV) the alpha reaches 8.9 fm from the gold centre
at b = 2 fm, inside the Woods-Saxon radius:

\snippet main_tutorial_forces_check.cpp nuclear

~~~{.txt}
b (fm)  Coulomb (deg)    +WoodsSaxon        +Yukawa  rmin (fm)
   2.0       129.3225        83.1825       129.2353       8.86
  30.0        16.0207        16.0207        16.0207      34.37
~~~

The angles are raw, at the stop (receding, 1000 fm); `rmin` is from the
Coulomb-only run. (The helper used here derives the event with a hand-written
loop `while(!ev.Inspect()) ev.DeriveDT(0.5);` to follow the smallest
distance; see the program.) The Woods-Saxon defaults are depth
\f$ V_0 \f$ = 50 MeV, radius \f$ R = 1.2\,(4^{1/3} + 197^{1/3}) \f$ fm =
8.886 fm (alpha + Au) and diffuseness 0.65 fm. The **Yukawa default strength
is only illustrative**, not fitted to anything: set it for your system. Both
act on **every** pair of the event, including alpha-alpha pairs in a
many-body event, with the same parameters. With a nuclear force present
`asymptote` is exact only when evaluated where the nuclear force is
negligible (it is, at 1000 fm); `asymptoteN` includes the nuclear part of
each pair (by quadrature of the orbit), also for a pair still approaching.

## Softening

An attractive Coulomb force diverges at \f$ r \to 0 \f$, and a head-on fall
into it fails. `coulombForce::SetSoftening(eps)` replaces
\f$ 1/r \f$ by \f$ 1/\sqrt{r^2+\varepsilon^2} \f$ (and the force by its exact
gradient), which keeps it finite; it changes the physics inside a few
\f$ \varepsilon \f$.

\snippet main_tutorial_forces_check.cpp soft

~~~{.txt}
VERLET softening 0 fm : ABORT     after  17550 steps, x =    -51.5 fm, |v| = 0.083650 c, |dE/E| = 1.0e-06
RK45   softening 0 fm : FAILED    after  17992 steps, x =     -nan fm, |v| = -nan c, |dE/E| = nan
RK45   softening 2 fm : INSPECTOR after  35992 steps, x =   1000.0 fm, |v| = 0.050000 c, |dE/E| = 1.9e-09
~~~

Without softening, the fixed-step `VERLET` steps over the singularity with a
large energy error, which the `"EDG"` abort catches. `RK45` instead shrinks
its substep inside a single `DeriveDT()` call, where no inspector is
evaluated; after `RK45MAXSUBSTEPS` (\f$ 10^6 \f$) substeps in one call it
gives up (with a message on stderr), `event::HasIntegratorFailed()` becomes
true for good, and `DeriveInspect()` ends with `kStopFailed` (`FAILED`). The
state is then not usable (here NaN). That takes about one second.

## Relativistic mode

`SetRelativistic()` integrates the momentum
\f$ p = \gamma m_0 v \f$ with \f$ dp/dt = F \f$. It works with `SEULER`,
`VERLET`, `RK4` and `RK45`; **`TAYLOR` and `oddevent` silently ignore it**.
The conserved energy is then `GetNetERel()`, not `GetNetE()`.

\snippet main_tutorial_forces_check.cpp rel

~~~{.txt}
non-relativistic  : angle 1.3991 deg, |v| 0.500000 c, |dE/E| Newton 4.4e-16, relativistic 3.1e-08
relativistic      : angle 1.2093 deg, |v| 0.500000 c, |dE/E| Newton 1.6e-08, relativistic 1.1e-14
~~~

The forces stay static (instantaneous Coulomb, no retardation, no magnetic
force). That is exact for a charge in the field of a fixed centre, as here.
For two fast **moving** charges the missing magnetic and retardation terms
are of the same order \f$ (v/c)^2 \f$ as the relativistic correction
itself, so the mode is not a consistent improvement there (\ref con_rel).
`asymptote` is non-relativistic and should not be used in this mode;
`asymptoteN` applies only its first-order correction there (flag
`RELATIVISTIC`), and its `kBrute` method integrates the relativistic copy.

## Your own force

Derive from \ref force and define `Evaluate(p1, p2)` (force on `p1` by `p2`,
amu c\f$^2\f$/fm) and `Potential(p1, p2)` (amu c\f$^2\f$). Two rules: the
force must be minus the gradient of the potential, or the energy checks are
meaningless; and `Evaluate(p2, p1)` must be `-Evaluate(p1, p2)`, because by
default each pair is evaluated once and the opposite force is given to the
other particle (otherwise call `event::SetPairSymmetric(kFALSE)`). The force
may depend on positions and on particle constants (`GetQ()`, `GetM0()`), not
on velocities (VERLET reuses the force of the previous step).

\snippet main_tutorial_forces_check.cpp myforce

Check both rules numerically before using it:

\snippet main_tutorial_forces_check.cpp myforceuse

~~~{.txt}
|F - (-grad U)| / |F| = 3.4e-11, |F12 + F21| / |F| = 0.0e+00
angle 38.0981 deg (unscreened: 51.01), |dE/E| = 2.0e-14
~~~

# 8. Recording tracks sparsely {#tut_record}

Program: `test/main_tutorial_record_check.cpp`.

A recording particle appends its state after every step, which for 100000
steps and several particles is a lot of memory. `SetRecordPeriod(n)` keeps
only the steps whose step count is a multiple of `n`:

\snippet main_tutorial_record_check.cpp record

~~~{.txt}
 period   steps   points       GetLastT        GetTime   |LastX - X| fm DCA err (fm)
      1   41378    41379        41378.0        41378.0           0.0000     4.97e-13
     10   41378     4138        41370.0        41378.0           0.4000     5.88e-11
    100   41378      414        41300.0        41378.0           3.8997     5.92e-09
~~~

- The initial state is **not** recorded unless you call `AddRecord()` before
  deriving (that is the "+1" in the number of points). `AddRecord()` always
  appends, whatever the period.
- With `n > 1`, `path::GetLastX()`, `GetLastV()`, `GetLastT()` are the last
  **recorded** step, not the final state, unless the step count happens to be
  a multiple of `n`. Use `particle::GetX()`, `GetV()`, `GetTime()` for the
  final state. (The project's own programs `main` and
  `test_rutherford_inspect` read the final velocity from `GetLastV()`, which
  is only right with period 1.)
- The parabolic refinement keeps the distance of closest approach accurate
  even from every 100th step of a smooth RK4 trajectory.

# Where next {#tut_next}

- \ref concepts "Concepts": the model behind each step of this tutorial.
- \ref validation "Validation": comparisons with known results.
- Class reference: \ref event, \ref particle, \ref force, \ref inspectorP,
  \ref inspectorE, \ref asymptote, \ref asymptoteN, \ref path, and `runinfo`.
