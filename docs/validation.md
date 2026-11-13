\page validation Physics validation

This page reports what the simulation has been checked against, and how well it
agrees. Every number below is copied from the output of the self-checking test
programs `test/main_validation_*_check.cpp`. They run under `ctest` with no
arguments, print their tables and return a non-zero exit code when a check fails.

| Test | Covers | Run time |
|---|---|---|
| `test_validation_rutherford_check` | finite-distance pitfall, Rutherford angle vs b for every integrator, brute-force cross-check, closest approach and `path::GetDCA`, relativistic mode | 0.8 s |
| `test_validation_recoil_check` | free (recoiling) target: lab angles, final speeds, momentum and energy | 0.4 s |
| `test_validation_crosssection_check` | Monte Carlo differential cross-section | 3.6 s |
| `test_validation_shortrange_check` | Coulomb + Woods-Saxon | 0.1 s |
| `test_validation_nbody_check` | three- and four-body conservation and convergence, a decoupled third body, convergence of the asymptotic state | 1.5 s |

Model background: \ref concepts. How to run: \ref tutorial.

# Two kinds of evidence {#val_kinds}

RooParticle is built for **three or more bodies**. For those there is no closed
form, so the checks are of two different kinds. Keep them apart when you quote a
number:

- **Exact two-body theory** (sections \ref val_pitfall to \ref val_shortrange).
  Two bodies with a pure \f$1/r\f$ interaction have closed forms: the Rutherford
  angle, the closest approach and the two-body kinematics. The relativistic
  fixed-centre deflection (Darwin) is also closed. These are the only *exact*
  references, and they test the integrators and the bookkeeping. They say
  nothing directly about many-body accuracy.
- **Numerical N-body evidence** (section \ref val_nbody): conservation laws,
  convergence order under step refinement, consistency with the two-body limit
  when a third body is decoupled, and convergence of the outgoing state with the
  stopping distance. These checks have no formula to agree with. They show that
  the numbers converge, and how fast.

Notation: \f$\kappa = k q_1 q_2\f$ with \f$k = 0.00154629\f$ amu c² fm / e²,
\f$\mu\f$ the reduced mass (the projectile mass for a fixed target), and
\f$a = \kappa/(\mu v_\infty^2)\f$.

# The finite-distance pitfall {#val_pitfall}

A run starts and stops at a finite distance, where the potential energy is not
negligible. At 1000 fm the Coulomb energy of alpha + Au is 4.9 % of the kinetic
energy. The textbook formula holds for the quantities **at infinity**. Those
follow exactly from the initial state through energy and angular momentum:

\f[ E = \tfrac12 \mu v_0^2 + \frac{\kappa}{r_0}, \qquad v_\infty = \sqrt{2E/\mu}, \qquad
    b_\infty = \frac{|\vec r_0 \times \vec v_0|}{v_\infty}, \qquad
    \tan\frac{\theta}{2} = \frac{\kappa}{\mu v_\infty^2 b_\infty} \f]

Alpha at 0.05c, nominal \f$b = 50\f$ fm, starting at \f$x_0 = -1000\f$ fm on a fixed Au:

| | nominal | asymptotic |
|---|---|---|
| impact parameter b (fm) | 50.000000 | 48.822869 |
| speed v (c) | 0.050000 | 0.051206 |
| Rutherford angle at these b, v (deg) | 52.082826 | 51.013865 |

- The Rutherford formula at the nominal values is off by **1.07 deg**.
- The raw angle between the start velocity and the velocity at the stop (1000 fm)
  is 50.945574 deg, off by 0.068 deg. This is smaller than the nominal error
  only because the tilt at the start and the residual deflection after the stop
  partly cancel.
- The angle between the asymptote of the initial state and that of the final
  state (class \ref asymptote) is 51.013865 deg. That is the asymptotic theory to
  1.0e-12 deg.

**Rule:** compare with theory using asymptotic quantities. Take them from
\ref asymptote for two bodies, or from energy and angular momentum computed by
hand. Never use the nominal \f$b\f$ and \f$v_0\f$. All theory values on this page
are computed by hand in the tests, like this:

\snippet main_validation_rutherford_check.cpp hand_theory

# Rutherford angle vs impact parameter {#val_rutherford}

Alpha on a fixed Au target, \f$x_0=-1000\f$ fm, stop when receding beyond 1000 fm.
There are 18 cases: \f$b\f$ from 0 to 500 fm, at \f$v_0 = 0.05c\f$ and \f$0.1c\f$.
The simulated angle is measured between the asymptotes of the initial and the
final state, and the theory is the hand formula above.

| Integrator | dt (fm/c) | max \f$\Delta\theta\f$ (deg) | max energy drift | tolerance (deg) |
|---|---|---|---|---|
| `TAYLOR` | 2 | 1.5e-1 | 9.0e-3 | 0.5 |
| `SEULER` | 2 | 4.4e-4 | 4.7e-6 | 2e-3 |
| `VERLET` | 2 | 4.1e-4 | 9.8e-12 | 2e-3 |
| `RK4` | 2 | 4.3e-9 | 4.0e-12 | 2e-8 |
| `RK45` | 20 (call step), default tol | 4.7e-8 | 1.0e-9 | 2e-7 |

Selected rows (simulation minus theory, deg):

| v0 | b_nom | b_inf | theta_R (deg) | TAYLOR | SEULER | VERLET | RK4 | RK45 |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 1 | 0.976 | 175.199243 | 2.1e-3 | 6.5e-6 | 8.8e-7 | -1.3e-13 | 1.3e-12 |
| 0.05 | 50 | 48.823 | 51.013865 | 7.9e-3 | 5.4e-5 | 2.0e-6 | 2.4e-13 | 6.6e-13 |
| 0.05 | 500 | 489.420 | 5.476601 | 1.7e-5 | -5.1e-6 | -1.6e-9 | 3.3e-14 | -3.8e-14 |
| 0.10 | 5 | 4.970 | 101.049986 | 1.5e-1 | 4.4e-4 | 4.1e-4 | -4.3e-9 | -4.7e-8 |
| 0.10 | 50 | 49.698 | 13.845758 | 2.5e-3 | 1.0e-5 | 1.8e-6 | -2.6e-12 | 6.8e-12 |

The worst case is always the hardest close encounter: 0.1c with b around 5 fm,
where the closest approach is about 12-15 fm. Head-on (\f$b=0\f$) gives exactly
180 deg for every integrator, by symmetry, so that row tests nothing about the
integrator.

**Independent of the \ref asymptote class.** For b = 5, 50, 500 fm the test also
integrates with RK45 (tolerance 1e-13) out to \f$10^9\f$ fm, both backward
(reversed velocity) and forward. It takes the angle between the two far
velocities directly:

| b_nom | theta_R (deg) | brute force (deg) | difference (rad) | (v - v_inf)/v |
|---|---|---|---|---|
| 5 | 156.324881000 | 156.324881000 | 2.7e-14 | -2.28e-8 |
| 50 | 51.013864754 | 51.013864754 | 1.3e-14 | -2.32e-8 |
| 500 | 5.476601115 | 5.476601115 | 5.2e-14 | -2.34e-8 |

Even at \f$10^9\f$ fm the speed is still below \f$v_\infty\f$ by
\f$\approx a/R = 2.3\times10^{-8}\f$. It agrees with
\f$v(R) = \sqrt{2(E-\kappa/R)/\mu}\f$ to 1.5e-13. Speeds converge like
\f$1/R\f$, directions faster (see \ref val_asym). `asymptote::GetThetaCM` of the
initial state equals the hand formula exactly (difference 0).

# Closest approach and path::GetDCA {#val_dca}

Closest approach for a repulsive \f$1/r\f$ potential:
\f$ r_{min} = a + \sqrt{a^2 + b_\infty^2} \f$, and head-on \f$r_{min} = 2a = \kappa/E\f$.
The comparison uses the full \f$E\f$, including \f$\kappa/r_0\f$. Using the nominal
kinetic energy instead gives 48.86 fm in place of 46.59 fm at 0.05c.

Alpha on fixed Au, dt = 1 fm/c. "traj" is the minimum of \f$r(t)\f$ of the run,
refined by a parabola through three steps. `GetDCA` is evaluated on the recorded
path (every step), and on every 25th sample of it (a 25 fm/c record spacing),
both discrete and interpolated. The bound is the largest possible discrete
sampling error, \f$\tfrac12 \ddot r\,(T/2)^2\f$.

| integ | v0 | b | d (fm) | traj - d | intp(1) - traj | disc(25) - traj | bound(25) | intp(25) - traj |
|---|---|---|---|---|---|---|---|---|
| VERLET | 0.05 | 0 | 46.58642 | -1.2e-6 | -7.1e-15 | 6.8e-4 | 2.2e-3 | -5.9e-8 |
| VERLET | 0.10 | 0 | 12.06827 | -1.8e-5 | -2.9e-12 | 2.2e-4 | 3.3e-2 | -1.2e-6 |
| VERLET | 0.05 | 20 | 53.69002 | -1.1e-6 | -9.9e-14 | 1.7e-3 | 2.2e-3 | -8.5e-8 |
| RK4 | 0.05 | 0 | 46.58642 | 1.4e-14 | -7.1e-15 | 6.8e-4 | 2.2e-3 | -5.9e-8 |
| RK4 | 0.10 | 0 | 12.06827 | 1.2e-11 | -3.0e-12 | 2.2e-4 | 3.3e-2 | -1.2e-6 |
| RK4 | 0.10 | 50 | 56.09697 | 1.5e-11 | -1.6e-11 | 3.8e-5 | 1.3e-2 | -2.0e-9 |

- The trajectory minimum agrees with the closed form to 1.7e-5 fm (VERLET) and
  1.5e-11 fm (RK4). That is the integration error.
- `GetDCA` with interpolation on the full record equals the trajectory minimum
  to 1.5e-11 fm.
- On the coarse record, the discrete DCA always lies inside its sampling bound
  (up to 1.7e-3 fm too large). The interpolated DCA is within 3.9e-5 of the bound
  (at most 1.2e-6 fm).
- The smallest DCA reported is 12.07 fm. Before the fix, `GetDCA` of a fixed
  point returned 0.

Free pairs, recorded on the same time grid, use the path-to-path
`GetDCA(path*)` (RK4, dt = 1), with \f$d = \kappa/E_{rel}\f$:

| pair | v0 | d (fm) | d nominal (fm) | traj - d | intp - traj | disc - traj |
|---|---|---|---|---|---|---|
| alpha + Au | 0.05 | 47.487421 | 49.854901 | -1.1e-13 | -1.4e-13 | 3.4e-7 |
| alpha + Au | 0.10 | 12.310293 | 12.463725 | 3.1e-10 | -4.5e-10 | 2.0e-5 |
| alpha + alpha | 0.02 | 15.227440 | 15.462900 | -1.1e-13 | -3.8e-13 | 6.4e-7 |

# Free (recoiling) target {#val_recoil}

Alpha at 0.05c from \f$x_0 = -1000\f$ fm on a target at rest at the origin, both
free. The textbook lab formulas assume the target is at rest **before the
collision, at** \f$t\to-\infty\f$. A target at rest at a finite distance was
already being pushed, so it had a velocity \f$v_{2,in}\f$ at \f$-\infty\f$:
\f$2.45\times10^{-5}c\f$ for Au and \f$3.09\times10^{-5}c\f$ for alpha
(6e-4 of \f$v_0\f$). The test finds \f$v_{2,in}\f$ by a backward RK45 integration
to \f$10^9\f$ fm and works in the frame S' where the target is at rest at
\f$-\infty\f$. In S', with \f$\theta\f$ the Rutherford angle at
\f$(b_\infty, v_\infty)\f$ of the relative motion:

\f[ \tan\theta_1 = \frac{\sin\theta}{\cos\theta + m_1/m_2},\qquad
    \theta_2 = \frac{\pi-\theta}{2},\qquad
    v_1' = v_\infty\frac{\sqrt{m_1^2+m_2^2+2m_1m_2\cos\theta}}{M},\qquad
    v_2' = \frac{2 m_1 v_\infty \sin(\theta/2)}{M} \f]

For equal masses, \f$\theta_1+\theta_2 = 90\f$ deg. The cases are Au (m = 197)
with b = 5, 20, 50, 100 fm, and alpha (m = 4) with b = 0.5, 1, 2, 5 fm. Each
integrator runs to a stop at 1000 fm, and the tail is extrapolated by
\ref asymptote. BRUTE is RK45 at tolerance 1e-13 run to \f$10^9\f$ fm, with no
asymptote class involved.

| integ | max \f$\Delta\theta\f$ lab (deg) | max \f$\Delta v/v\f$ | max \f$\Delta P/P\f$ | max \f$\Delta E/E\f$ | equal masses: \f$\theta_1+\theta_2-90\f$ (deg) |
|---|---|---|---|---|---|
| VERLET dt = 1 | 3.1e-4 | 1.2e-5 | 1.7e-14 | 2.5e-13 | 6.4e-13 |
| RK4 dt = 2 | 7.9e-8 | 2.9e-9 | 2.5e-14 | 1.4e-10 | 1.1e-8 |
| RK45 default tol | 1.5e-7 | 3.2e-9 | 6.1e-15 | 1.4e-9 | 8.9e-8 |
| BRUTE | 2.3e-12 | 4.8e-14 | 1.3e-14 | 2.2e-13 | 7.0e-13 |

Example values: Au, b = 20 fm gives \f$\theta_{CM}\f$ = 101.1540, \f$\theta_1\f$ =
100.0083, \f$\theta_2\f$ = 39.4230 deg. Alpha, b = 1 fm gives \f$\theta_{CM}\f$ =
102.0276 and \f$\theta_1 = \theta_{CM}/2\f$ = 51.0138 deg. The worst cases are
alpha + alpha at b of 0.5-1 fm (closest approach 2.6 fm). Total momentum is
conserved to round-off (pair forces are equal and opposite). Taking the textbook
lab angle at the nominal \f$(v_0, b)\f$ and comparing it with the raw velocity at
the stop is wrong by up to **1.4 deg**.

# Relativistic mode {#val_rel}

With `event::SetRelativistic(kTRUE)` the integrators solve \f$d\vec p/dt = \vec F\f$
with \f$H = \sqrt{p^2+m^2} + \kappa/r\f$ (c = 1). For a fixed centre this has an
exact orbit, \f$u'' + \Gamma^2 u = -\kappa W/L^2\f$, with \f$u = 1/r\f$,
\f$\Gamma^2 = 1-\kappa^2/L^2\f$, \f$W\f$ the total energy and \f$L = |\vec r\times\vec p|\f$.
It gives the deflection (Darwin):

\f[ \theta = \pi - \frac{2}{\Gamma}\arctan\frac{\Gamma\, p_\infty L}{\kappa W},
    \qquad p_\infty = \sqrt{W^2-m^2} \f]

\snippet main_validation_rutherford_check.cpp darwin

The simulation runs RK45 at tolerance 1e-13, backward and forward to
\f$10^9\f$ fm or more. "NR" is \f$2\arctan(\kappa/(m_0 v^2 b))\f$ and "pv" is
\f$2\arctan(\kappa/(p v b))\f$, both evaluated at the asymptotic \f$v\f$, \f$p\f$ and \f$b\f$.

| projectile | v0 | b (fm) | kappa/L | Darwin (deg) | sim - Darwin (rad) | NR - sim (deg) | pv - sim (deg) |
|---|---|---|---|---|---|---|---|
| alpha | 0.01 | 1057.857 | 0.006 | 58.617886 | 2.2e-13 | 3.8e-3 | 1.2e-3 |
| alpha | 0.05 | 42.263 | 0.029 | 58.592228 | -5.5e-14 | 9.5e-2 | 3.0e-2 |
| alpha | 0.05 | 6.538 | 0.187 | 148.583149 | -2.4e-14 | 6.8e-2 | 2.7e-2 |
| alpha | 0.30 | 1.121 | 0.173 | 57.643624 | -1.9e-14 | 3.57 | 1.10 |
| positron | 0.50 | 1336.070 | 0.289 | 55.839752 | -5.3e-13 | 10.8 | 3.13 |
| positron | 0.90 | 679.596 | 0.159 | 18.167896 | -7.9e-14 | 26.1 | 1.80 |
| positron | 0.90 | 207.554 | 0.520 | 48.699925 | -1.8e-14 | 58.2 | 11.0 |
| positron | 0.99 | 55.513 | 0.572 | 46.225992 | -2.1e-14 | 107 | 13.7 |

- The relativistic simulation reproduces the exact relativistic angle to 5e-13
  rad at every speed. The relativistic energy (`event::GetNetERel`) is conserved
  to 6.5e-13.
- At low speed it reproduces the non-relativistic Rutherford angle up to the
  expected \f$O(v^2)\f$ correction. The difference is below \f$0.61\,\theta v^2\f$
  (0.0038 deg at 0.01c, 0.095 deg at 0.05c).
- At high speed (information) the non-relativistic formula moves away fast: 3.6
  deg at 0.3c, and tens of degrees for a positron at 0.9c. The "pv" form is
  better but still off by up to 14 deg when \f$\kappa/L\f$ is large.
- Retardation and radiation are not in the model. These checks validate the
  integration of the stated Hamiltonian, not relativistic electrodynamics.

# Differential cross-section {#val_xsec}

Rutherford: \f$\frac{d\sigma}{d\Omega} = \left(\frac{a}{2}\right)^2 \frac{1}{\sin^4(\theta/2)}\f$.
A beam uniform over the disc \f$b_\infty < b_{max}\f$ is **uniform in \f$b^2\f$, not in
\f$b\f$**. The expected number of events in a bin \f$[\theta_1,\theta_2]\f$ is the exact
integral \f$N\,(b(\theta_1)^2 - b(\theta_2)^2)/b_{max}^2\f$ with
\f$b(\theta) = a\cot(\theta/2)\f$, so there is no bin-centre approximation.

Set-up: alpha at 0.05c on fixed Au, 20000 events, fixed seed. Every event starts
1000 fm from the target, so every event has the same energy and
\f$b_\infty = b\,v_0/v_\infty\f$ exactly. Then \f$a\f$ = 23.2932 fm and
\f$b_{max}\f$ = 147.068 fm (\f$\theta_{min}\f$ = 18 deg), and the beam area is
679.5 b. Each event runs RK45 with the default tolerance and a 400 fm/c call
step, and stops with the inspectors `"VRG"` > 0 and `"DSTG"` > 1000 (about 111
RK45 substeps per event).

| theta (deg) | N obs | N exp | pull | measured (mb/sr) | bin average (mb/sr) | ratio |
|---|---|---|---|---|---|---|
| 20-25 | 5898 | 5928.7 | -0.40 | 9.553e5 | 9.603e5 | 0.995 |
| 25-30 | 3238 | 3220.1 | 0.31 | 4.346e5 | 4.322e5 | 1.006 |
| 30-40 | 3242 | 3200.7 | 0.73 | 1.753e5 | 1.731e5 | 1.013 |
| 40-50 | 1487 | 1479.9 | 0.18 | 6.523e4 | 6.492e4 | 1.005 |
| 50-60 | 753 | 802.2 | -1.74 | 2.852e4 | 3.038e4 | 0.939 |
| 60-75 | 655 | 653.0 | 0.08 | 1.468e4 | 1.464e4 | 1.003 |
| 75-90 | 353 | 350.4 | 0.14 | 7375 | 7320 | 1.007 |
| 90-110 | 288 | 255.7 | 2.02 | 4553 | 4043 | 1.126 |
| 110-135 | 137 | 159.9 | -1.81 | 2029 | 2368 | 0.857 |
| 135-180 | 100 | 86.1 | 1.50 | 1846 | 1589 | 1.162 |

- Uniform in \f$b^2\f$: \f$\chi^2\f$ = 13.47 for 10 bins, p = 0.20, largest pull 2.0.
  The deviations are statistical: every event agrees with the Rutherford angle
  at its own \f$b_\infty\f$ to 4.5e-9 rad.
- Uniform in \f$b\f$ (the wrong weighting, 2000 events): \f$\chi^2\f$ = 2804 for 10
  bins. The test is sensitive to the weighting.

# Short-range force: Coulomb + Woods-Saxon {#val_shortrange}

Alpha at 0.12c on fixed Au (\f$E_\infty\f$ = 27.06 MeV, \f$a\f$ = 4.206 fm), using
`coulombForce` alone and `coulombForce` + `woodsSaxonForce` with the defaults
(V0 = 50 MeV, R = 8.886 fm, diffuseness 0.65 fm). At 1000 fm the nuclear term is
\f$e^{-1525}\f$, so the asymptotes at the start and at the stop are pure Coulomb.
There is no closed form with the nuclear force.

| b_nom (fm) | Coulomb r_min (fm) | theta_R (deg) | C+WS (deg) | C+WS - theta_R (deg) |
|---|---|---|---|---|
| 1 | 8.53 | 153.3598 | 126.9625 | -26.4 |
| 3 | 9.36 | 109.2291 | 48.6035 | -60.6 |
| 5 | 10.72 | 80.3780 | 14.0039 | -66.4 |
| 6 | 11.51 | 70.2870 | 92.2863 | +22.0 |
| 7 | 12.35 | 62.2121 | 58.5888 | -3.62 |
| 9 | 14.11 | 50.2811 | 50.0607 | -0.220 |
| 10 | 15.02 | 45.7953 | 45.7395 | -5.6e-2 |
| 12 | 16.87 | 38.7814 | 38.7780 | -3.4e-3 |
| 15 | 19.72 | 31.4522 | 31.4521 | -4.6e-5 |
| 20 | 24.56 | 23.8495 | 23.8495 | -3.0e-8 |
| 30 | 34.37 | 16.0279 | 16.0279 | -3.5e-10 |

- Coulomb alone matches Rutherford at every b to 9.7e-11 rad.
- With the nuclear force, scattering at large b stays Rutherford: below 6.2e-12
  rad for b of 30 fm and up.
- Small b is changed strongly. The angle falls far below Rutherford, then jumps
  above it at b = 6 fm: the orbit is refracted by the attractive well, the
  nuclear-rainbow behaviour.
- (Information) The deviation exceeds 1e-3 of the angle up to b = 10 fm, where
  the Coulomb closest approach is 15.0 fm, i.e. R + 9.4 diffuseness. It exceeds
  0.1 deg up to b = 9 fm (14.1 fm, R + 8.0 diffuseness). It falls off like the
  Woods-Saxon tail \f$e^{-(r-R)/a}\f$, by about x4 per fm of b.
- Energy is conserved to 2.5e-10 with RK45 at tolerance 1e-10, over the whole
  run with both forces.
- With no formula available, the combined-force angle is cross-checked by two
  integrators. RK45 at tolerance 1e-10 vs 1e-12 differs by up to 1.6e-8 deg.
  VERLET at dt = 0.2 fm/c differs from RK45 by up to 2.7e-4 deg, and the error
  ratio between dt = 0.4 and 0.2 is 4.00 (second order). VERLET's energy error
  inside the well reaches 2.9e-5 at dt = 0.2.

# N-body: numerical evidence {#val_nbody}

## Conservation and order of convergence

There are three configurations: all bodies free, a fixed target with two
projectiles that also repel each other, and four free bodies:

\snippet main_validation_nbody_check.cpp configs

The smallest pair distance during the runs is 20.5 fm (C1), 55.7 fm (C2) and
22.3 fm (C3). The columns are:

- dE: the largest \f$|E-E_0|/|E_0|\f$ over the run.
- dP and dL: the final drift divided by \f$\sum|\vec p_i|\f$ and \f$\sum|\vec x_i\times\vec p_i|\f$.
- dx: the final position error against a reference, RK4 at dt = 0.25. RK45 at
  tolerance 1e-13 agrees with that reference to 2.9e-10 fm.
- p: the order from step halving.

In C2 momentum is not conserved (the fixed target absorbs it), but angular
momentum about the target is.

| config | integ | dt | dE | dP | dL | dx (fm) | p |
|---|---|---|---|---|---|---|---|
| C1 3 free | TAYLOR | 1 | 4.0e-5 | 3.4e-15 | 5.5e-5 | 3.3e-3 | 1.00 |
| C1 3 free | SEULER | 1 | 3.7e-5 | 1.4e-15 | 3.1e-15 | 2.1e-4 | 1.00 |
| C1 3 free | VERLET | 1 | 4.1e-8 | 1.9e-15 | 3.7e-15 | 1.6e-6 | 2.00 |
| C1 3 free | RK4 | 4 | 3.8e-12 | 1.4e-15 | 4.2e-15 | 4.2e-10 | 3.99 |
| C2 Au + 2 alpha | TAYLOR | 1 | 2.8e-4 | - | 1.0e-4 | 2.1e-1 | 1.00 |
| C2 Au + 2 alpha | SEULER | 1 | 1.1e-4 | - | 4.7e-14 | 1.4e-3 | 1.01 |
| C2 Au + 2 alpha | VERLET | 1 | 2.0e-8 | - | 8.1e-14 | 2.6e-5 | 2.00 |
| C2 Au + 2 alpha | RK4 | 4 | 9.7e-14 | - | 1.2e-14 | 3.1e-10 | 3.93 |
| C3 4 free | TAYLOR | 1 | 1.5e-4 | 6.2e-15 | 1.8e-4 | 1.4e-1 | 1.00 |
| C3 4 free | SEULER | 1 | 8.0e-5 | 6.1e-15 | 6.2e-14 | 2.0e-3 | 1.00 |
| C3 4 free | VERLET | 1 | 1.5e-8 | 1.4e-14 | 4.6e-14 | 2.2e-5 | 2.00 |
| C3 4 free | RK4 | 4 | 1.8e-13 | 3.0e-15 | 3.1e-14 | 1.6e-9 | 4.04 |

RK45 (call step 200 fm/c, so that the tolerance, not the call step, sets the substep):

| config | tolerance | dE | dL | dx (fm) | substeps |
|---|---|---|---|---|---|
| C1 | 1e-6 / 1e-8 / 1e-10 | 2.4e-7 / 1.6e-9 / 1.7e-11 | 2.6e-9 / 2.8e-11 / 8.9e-13 | 7.2e-5 / 6.8e-7 / 5.5e-9 | 103 / 116 / 162 |
| C2 | 1e-6 / 1e-8 / 1e-10 | 1.1e-8 / 9.9e-10 / 3.5e-11 | 1.2e-10 / 1.1e-10 / 2.2e-12 | 6.2e-6 / 1.8e-6 / 3.2e-8 | 200 / 217 / 271 |
| C3 | 1e-6 / 1e-8 / 1e-10 | 1.3e-8 / 1.0e-9 / 2.2e-11 | 5.1e-9 / 1.9e-10 / 2.4e-12 | 4.8e-5 / 3.2e-6 / 7.5e-8 | 200 / 220 / 294 |

- The measured orders are 1, 1, 2 and 4 for TAYLOR, SEULER, VERLET and RK4 in
  all three configurations. The check range is \f$\pm0.3\f$ and the measured
  values lie within 0.0-0.07 of the expected ones.
- Momentum is conserved to round-off (below 1e-13) by every integrator, because
  pair forces are equal and opposite.
- SEULER and VERLET conserve angular momentum to round-off (symplectic, central
  pair forces). RK4 and RK45 conserve it to their accuracy. TAYLOR drifts at
  first order (1e-4).
- In RK45, `DeriveDT(dt)` caps the substep at `dt`. With a small call step
  (4 fm/c) the results were identical at every tolerance from 1e-6 to 1e-13.
  Give RK45 a call step larger than the substep it needs.

## A far, nearly decoupled third body {#val_decouple}

Alpha (b = 50 fm, 0.05c) on fixed Au, plus a proton at rest at \f$(0, D, 0)\f$.
The pair must reproduce the exact two-body Rutherford angle to within the size
of the perturbation. At the pair, the proton's field is nearly uniform,
\f$F = k q_\alpha q_p/D^2\f$, over the run time \f$T\f$ = 41380 fm/c. To first
order it changes the outgoing direction in two ways:

- directly, by the impulse \f$F T/(m v_0)\f$;
- through its torque \f$x F\f$ about the target, which changes \f$L\f$ and so \f$b\f$.

The bound is

\f[ \delta\theta_{est} = \frac{F}{m v_0}\left[T + \left|\frac{d\theta}{db}\right|\int |x(t)|\,dt\right] \f]

with \f$|d\theta/db|\f$ = 0.01592 rad/fm and \f$\int|x|dt\f$ = 1.588e7 fm²/c. The
torque term is 6.1 times the direct one. Without it, the deviation would be 3
times the naive estimate.

| D (fm) | estimate (rad) | direction change vs 2-body (rad) | theta - theta_R (rad) | ratio | slope |
|---|---|---|---|---|---|
| 1e4 | 4.550e-5 | 1.942e-5 | 1.942e-5 | 0.427 | |
| 3e4 | 5.056e-6 | 2.161e-6 | 2.161e-6 | 0.427 | 1.999 |
| 1e5 | 4.550e-7 | 1.941e-7 | 1.941e-7 | 0.427 | 2.001 |
| 3e5 | 5.056e-8 | 2.156e-8 | 2.156e-8 | 0.426 | 2.001 |
| 1e6 | 4.550e-9 | 1.940e-9 | 1.939e-9 | 0.426 | 2.000 |

The three-body result reproduces Rutherford within the estimate at every \f$D\f$
(ratio 0.43, the same at every \f$D\f$, so the estimate is tight and not merely a
bound). The deviation scales as \f$D^{-2.000}\f$. Without the proton, the
two-body run gives Rutherford to 6.7e-15 rad.

## Asymptotic state of three or more bodies {#val_asym}

For three or more bodies no formula gives the state at infinity. The only way
is to integrate far. This test uses RK45 at tolerance 1e-13 out to
\f$10^{12}\f$ fm, with a call step proportional to the smallest pair distance:

\snippet main_validation_nbody_check.cpp far_loop

\f$\delta(R)\f$ is the largest angle, over the bodies, between the velocity
direction at stopping distance \f$R\f$ (the smallest pair distance) and the
direction at \f$R = 10^{12}\f$ fm. The fit is \f$\delta = C/R^k\f$ over
\f$10^4\f$ to \f$10^8\f$ fm.

| system | 1e3 | 1e4 | 1e5 | 1e6 | 1e7 | 1e8 | k | C | R for 1e-6 rad (fm) |
|---|---|---|---|---|---|---|---|---|---|
| alpha + fixed Au | 5.9e-4 | 5.7e-6 | 5.6e-8 | 5.6e-10 | 5.7e-12 | 5.8e-14 | 2.00 | 5.6e2 | 2.4e4 |
| alpha + free Au, relative velocity | 6.1e-4 | 5.7e-6 | 5.8e-8 | 5.7e-10 | 5.7e-12 | 6.1e-14 | 2.00 | 5.6e2 | 2.4e4 |
| alpha + free Au, lab | 2.5e-2 | 2.4e-3 | 2.4e-4 | 2.4e-5 | 2.4e-6 | 2.4e-7 | 1.00 | 24 | 2.4e7 |
| C1 3 free, lab | 1.3e-3 | 1.3e-4 | 1.3e-5 | 1.3e-6 | 1.3e-7 | 1.3e-8 | 1.00 | 1.3 | 1.3e6 |
| C2 fixed Au + 2 alpha, lab | 5.9e-4 | 2.9e-5 | 2.6e-6 | 2.6e-7 | 2.6e-8 | 2.6e-9 | 1.01 | 0.31 | 2.7e5 |
| C3 4 free, lab | 5.2e-3 | 5.4e-4 | 5.4e-5 | 5.4e-6 | 5.4e-7 | 5.4e-8 | 1.00 | 5.4 | 5.4e6 |

What this tells a user:

- With **two bodies**, the relative direction converges as \f$1/R^2\f$.
  Physically, the force is radial and nearly parallel to the velocity, so only
  the offset \f$b/R\f$ deflects it. The \ref asymptote class makes this moot:
  it is exact.
- With **three or more bodies**, the lab directions converge only as
  **\f$1/R\f$**. The force between two receding bodies is along their
  separation, which is in general not along either velocity. This also holds
  for two free bodies seen in the lab, where the relative speed still changes
  like \f$1/R\f$. Speeds converge as \f$1/R\f$ in every case (the remaining
  potential energy).
- So the stopping distance must be \f$R \gtrsim C/\delta\f$ for an accuracy
  \f$\delta\f$. For these systems that means \f$10^5\f$ to \f$10^7\f$ fm for
  1e-6 rad, and one more decade per extra digit. The usual 1000 fm stop leaves
  errors of 1e-3 to 5e-3 rad (0.07-0.3 deg) in the N-body outgoing directions.
  Integrating that far is cheap with RK45 and a call step proportional to the
  distance: about 2% growth per call, so a few hundred calls per decade.

# Limits {#val_limits}

- The closed forms used are **two-body** and, except for the Darwin angle,
  **non-relativistic**. The N-body results are numerical evidence (convergence,
  conservation, consistency), not comparisons with exact theory.
- The tolerances in the tests are about 3-7 times the measured values **for these
  specific set-ups** (alpha on Au at 0.05-0.12c, the given steps). They are not
  general accuracy guarantees. Closer encounters, higher speeds or larger steps
  give larger errors. Measure them for your own system the same way.
- The relativistic mode is checked as mechanics (\f$d\vec p/dt = \vec F\f$ with a
  static field). Retardation and radiation are not modelled.
- `woodsSaxonForce` is checked for consistency (energy conservation, integrator
  agreement, a Rutherford limit at large b), not against nuclear data.
- `oddevent` is not covered here.
