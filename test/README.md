# Test programs

Every `test/main_<name>.cpp` is built as the executable `test_<name>` in the
build directory. Files listed in `RP_TEST_EXCLUDE` in `test/CMakeLists.txt`
are skipped (at present `main_rutherford_track.cpp`, which does not compile).

## Self-checking tests (run by ctest)

Every `main_*_check.cpp` and `main_integrator_compare.cpp` takes no
arguments, prints a table of what it checked and exits with a non-zero code
when a check fails. CMake registers them in `ctest` automatically (re-run
`cmake` after adding one); CI runs `ctest` on every push.

| Kind | Programs | What they check |
|---|---|---|
| Library | `asymptote_check`, `asymptoteN_check`, `event_check`, `force_check`, `inspector_check`, `path_check`, `integrator_compare` | the two-body and N-body asymptotes, integrators (order, RK45, relativistic mode, conservation, pair symmetry), ownership and copying, forces against their potentials, every inspector method, closest approach, terminating conditions |
| Tutorial | `tutorial_*_check` (first, stop, physics, multibody, asymptoteN, montecarlo, accuracy, forces, record) | the programs shown in the Tutorial page of the manual; each checks the results the tutorial states |
| Validation | `validation_*_check` (rutherford, recoil, crosssection, shortrange, nbody) | physics against exact two-body results and numerical evidence for three and more bodies, reported on the Validation page of the manual |

`test_tutorial_montecarlo_check` also accepts `<title> <jobIndex>`, the form
used by `macro/runmulti.py`, and writes `Data/<title>_<jobIndex>.root`.

## Simulation programs (not in ctest)

Longer runs with command-line options, some of which open ROOT canvases:
`rutherford`, `rutherford_inspect` (`-h` lists its options), `rutherford_inspect_odd`,
`rutherford_inspect_oneparticle_showtrack`, `rutherford_inspect_oneparticle_tracksave`,
`rutherford_recorder`, `3body_oneline`, `3body_showtrack`, `3body_uniform`,
`collider_HeHe` (`[title] [jobIndex [seed]]`), `colliding_inspect`,
`colliding_inspect_odd`, `preEventTest`. They write to `Data/` in the working
directory. They are examples of use, not tests: they do not check their
results.
