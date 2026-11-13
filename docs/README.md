# How to write docs here

Not part of the generated manual (excluded in `Doxyfile.in`).
Build: `cmake <src> && make doc`, open `<build>/docs/html/index.html`.
Warnings: `<build>/docs/doxygen-warnings.log`. Formulas use MathJax from a CDN
(needs internet when viewing).

## Pages

`mainpage.md`, `tutorial.md`, `concepts.md`, `validation.md` are Markdown pages.
Each starts with `\page <id> <Title>` (the id is stable; never rename it), the
main page with `\mainpage`. They show up under "Related Pages".
A new page = a new `docs/<name>.md` with its own `\page` line, linked from `mainpage.md`.

## Links

- Page: `\ref tutorial "Tutorial"`
- Class: `\ref event`, `\ref inspectorP`; member: `\ref event::AddForce`
- Class names in running text are linked automatically.

## Formulas

Inline `\f$ E = \frac{1}{2} m v^2 \f$`, display:

    \f[ \tan\frac{\theta}{2} = \frac{Z_1 Z_2 e^2}{2 E b} \f]

## Compilable examples (`EXAMPLE_PATH = test`)

Whole file: `\include main_rutherford.cpp`. Part of a file: mark the region in the
`.cpp` with `//! [setup]` ... `//! [setup]` and write `\snippet main_rutherford.cpp setup`.
The examples are the test programs, so they stay compilable.

## Comment convention for headers

Use `///` with `\brief`, `\param`, `\return`. Put the comment on the lines before
the declaration, in the header. Use `///<` for a member documented on its own line.
State units in the text (fm, amu, c, e, fm/c).

```cpp
/// \brief Point-like charge interacting through the Coulomb force.
///
/// Mass in amu, charge in e, position in fm, velocity in c.
class EMparticle : public particle {
public:
  /// \brief Advance the particle by one time step.
  /// \param dt time step in fm/c
  /// \return false if the particle did not move (invincible)
  bool Step(double dt);

  double mass; ///< mass in amu
};
```
