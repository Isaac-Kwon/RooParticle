# RooParticle

Classical scattering of three or more charged bodies, built on ROOT (CERN).

[![CI](https://github.com/Isaac-Kwon/RooParticle/actions/workflows/ci.yml/badge.svg)](https://github.com/Isaac-Kwon/RooParticle/actions/workflows/ci.yml)

RooParticle is a C++ library and a set of programs for many-body scattering
events. Particles (point masses, usually charged nuclei) interact through
pairwise forces: Coulomb, nuclear potentials, or forces you write yourself.
An event is integrated in time until terminating conditions that you choose
are met. The run must stop at a finite distance, so asymptotic velocities are
computed with an error estimate. Results and run metadata are written to a
ROOT `TTree`. Two-body scattering (Rutherford) is the special case with an
exact answer, and is used for validation.

The manual is generated with Doxygen (`make doc`, see below). It has a
tutorial in which every code block is a tested program, a concepts page and
a validation page. Open work is listed in [TODO.md](TODO.md).

## Units

| Quantity | Unit |
|---|---|
| length | fm |
| mass | amu |
| velocity | c (speed of light) |
| charge | e |
| time | fm/c |
| energy | amu c² (= 931.494 MeV) |

## Build, test, documentation

Requires ROOT, CMake (3.16 or newer) and a C++ compiler (Linux or WSL).
Doxygen is optional and only needed for the manual.

```bash
source <ROOT directory>/bin/thisroot.sh
mkdir build && cd build
cmake <RooParticle directory>
make -j8
ctest --output-on-failure    # self-checking tests, tutorial and validation programs
make doc                     # manual in build/docs/html/index.html
```

The library `librooparticle.a`, the program `main` and one executable
`test_<name>` per `test/main_<name>.cpp` are built in the build directory.
Programs write their ROOT files to `Data/` in the working directory.

## A minimal example

An alpha particle on a fixed gold-like nucleus, integrated with velocity
Verlet for 40000 steps of 1 fm/c. The setup is copied verbatim from
`test/main_tutorial_first_check.cpp`,
which `ctest` compiles and runs; the full program also checks the angle and
the energy conservation.

```cpp
#include "event.hpp"
#include "particle.hpp"
#include "force.hpp"

// The force is not owned by the event: keep it alive as long as the event.
coulombForce coulomb;

event ev;
ev.AddForce(&coulomb);

// EMparticle(mass [amu], charge [e], position [fm], velocity [c], invincible, record)
// The event owns the particles: create them with new and never delete them yourself.
particle* target = ev.AddParticle(
    new EMparticle(197., 79., TVector3(0., 0., 0.), TVector3(0., 0., 0.), true, false));
particle* alpha  = ev.AddParticle(
    new EMparticle(4., 2., TVector3(-1000., 50., 0.), TVector3(0.05, 0., 0.), false, false));

ev.SetIntegrator("VERLET");
ev.DeriveDTN(40000, 1.);             // 40000 steps of dt = 1 fm/c

TVector3 v = alpha->GetV();          // final velocity [c]
```

An *invincible* particle is fixed: it acts on the others but never moves. To
stop automatically instead of after a fixed number of steps, add
*inspectors* (terminating conditions) and call `DeriveInspect()`.

## Features

**Integrators** (`event::SetIntegrator()`):

| Name | Scheme | Order | Force evaluations per step |
|---|---|---|---|
| `TAYLOR` | original scheme (default) | 1 | 1 |
| `SEULER` | symplectic Euler | 1 | 1 |
| `VERLET` | velocity Verlet, symplectic | 2 | 1 |
| `RK4` | classical Runge-Kutta | 4 | 4 |
| `RK45` | Dormand-Prince 5(4), adaptive substeps with tolerances | 5 | 6 per substep |

A relativistic mode (`SetRelativistic()`, momentum integration with static
forces) works with all of them except `TAYLOR`. `oddevent` is a variant that
moves each particle right after its own force is computed.

**Forces** (`force.hpp`): `coulombForce` (with optional softening),
`universalForce` (gravity), `yukawaForce`, `woodsSaxonForce`, and your own
subclass of `force`.

**Terminating conditions**: pair inspectors (`inspectorP`: angle between
relative position and velocity, distance, radial velocity, speed, step count)
and event inspectors (`inspectorE`: energies, energy drift, time, smallest
distance, largest acceleration, leaving a volume, and `ASYD`/`ASYV`, the
remaining deflection or velocity change from now to infinity). They are
combined with AND (default) or OR; abort inspectors and a step cap end an
event early, and `event::GetStopReason()` tells why it ended.

**Asymptotes**: `asymptote` gives the exact incoming and outgoing state of a
two-body pure 1/r collision from any instant. `asymptoteN` gives the
asymptotic velocities of any number of bodies, with error estimates, bound
clusters and reliability flags.

**Many events**: template events and their copies, run metadata (seed,
integrator, inspectors, git revision) in the ROOT file through `runinfo`,
the parallel job runner `macro/runmulti.py` and the plot macro
`macro/plot_rutherford.C`.

**Tests**: self-checking programs for the library, the tutorial and the
physics validation, run by `ctest` and in CI; see
[test/README.md](test/README.md).

## Documentation

After `make doc`, open `build/docs/html/index.html`:

- *Tutorial*: from the first event to stopping conditions, two-body
  asymptotes, events with three and more bodies and their asymptotic
  velocities, Monte Carlo runs with parallel jobs, accuracy, other forces and
  track recording.
- *Concepts*: units and constants, ownership and copying of events, the
  integrators, every terminating condition, the finite-distance problem and
  the asymptote methods, forces, reproducibility.
- *Validation*: comparisons with exact two-body results and numerical
  evidence for three and more bodies.
- The class reference, generated from the headers in `include/`.

The sources of these pages are in `docs/`; `docs/README.md` explains how to
write them.
