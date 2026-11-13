\mainpage RooParticle

RooParticle is a classical scattering simulator built on ROOT for three or
more bodies. Charged particles interact through pairwise forces (Coulomb,
optionally nuclear or your own), are integrated in time until terminating
conditions that you choose are met, and the results (velocities, angles,
closest approaches, with the run parameters) go to a ROOT `TTree`. Because a
run must stop at a finite distance, \ref asymptoteN corrects the final
velocities of any number of bodies to their values at infinity, with an
error estimate (\ref asymptote does it exactly for two bodies).
Units: length fm, mass amu, velocity c, charge e, time fm/c, energy
amu c\f$^2\f$ (931.494 MeV).

## A minimal example

An alpha particle on a fixed gold-like nucleus (from
`test/main_tutorial_first_check.cpp`, which `ctest` compiles and runs):

\snippet main_tutorial_first_check.cpp setup

\snippet main_tutorial_first_check.cpp derive

\snippet main_tutorial_first_check.cpp result

## Where to go

- \ref tutorial "Tutorial": from this first event to stopping conditions,
  exact two-body asymptotes, events with three and more bodies and their
  asymptotic velocities, Monte Carlo runs with parallel jobs, accuracy, other
  forces and track recording. Every code block is a tested program.
- \ref concepts "Concepts": units and constants, ownership and copying,
  integrators, terminating conditions in full, the finite-distance problem,
  forces, reproducibility.
- \ref validation "Validation": physics checks against known results.
- Class reference: \ref event, \ref particle, \ref force, \ref inspectorP,
  \ref inspectorE, \ref asymptote, \ref asymptoteN (full list under "Classes").
- Build: `cmake <source> && make && ctest`; this manual: `make doc`
  (\ref tut_build). Writers of these pages: `docs/README.md` in the repository.
