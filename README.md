<div align="center">

# marl

Explicit ODE solvers in plain Julia: forward Euler, classic Runge-Kutta 4 and an adaptive Dormand-Prince 5(4) pair. No dependencies.

![License](https://img.shields.io/badge/LICENSE-MIT-blue?style=for-the-badge)  ![Language](https://img.shields.io/badge/JULIA-9558b2?style=for-the-badge&logo=julia&logoColor=white)  ![Tests](https://img.shields.io/badge/TESTS-21-brightgreen?style=for-the-badge)

Refer to [LICENSE](LICENSE) for the legal bits.

</div>

## What's in it

- `euler` and `rk4` take a fixed step `dt`, and the last step is shortened so the solution lands exactly on the end of the interval.
- `rk45` is Dormand-Prince 5(4). It picks its own step from absolute and relative tolerances (`abstol`, `reltol`), reuses the last stage of each accepted step as the first stage of the next, and gives up with an error after `maxsteps` attempts.
- Everything is written with broadcasting, so the same code handles scalars and vectors.
- The returned `Solution` holds the time points and states, counts accepted and rejected steps, and can be called like a function to interpolate linearly between points. Asking outside the solved interval throws a `DomainError`.

The right-hand side is always `f(t, u)` and returns the derivative. It isn't mutated in place.

```julia
using Marl

oscillator(t, u) = [u[2], -u[1]]

sol = rk45(oscillator, [1.0, 0.0], (0.0, 10.0); abstol = 1e-10, reltol = 1e-10)

last(sol)
sol(2.5)
sol.accepted, sol.rejected
```

## 🏗️ Architecture

```
 src/Marl.jl  (module Marl, exports Solution euler rk4 rk45)
 ┌───────────────────────────────────────────────────────────┐
 │ Solution{T}  t::Vector{Float64}  u::Vector{T}             │
 │              accepted  rejected     (sol)(t) -> lerp      │
 │                                                           │
 │ fixed(step, f, u0, tspan, dt)  <──  euler_step            │
 │   loop with a shortened final step   rk4_step             │
 │                                                           │
 │ rk45:  tableau A, C, B4 (const tuples)                    │
 │   stages 1..6 -> fifth / fourth order estimates           │
 │   error_norm -> accept or reject -> rescale h             │
 └───────────────────────────────────────────────────────────┘
 test/runtests.jl      Project.toml
```

- `euler` and `rk4` share one driver, `fixed`, and differ only in the single-step function they pass.
- `rk45` rescales the step by `0.9 * err^(-1/5)`, clamped between 0.2x and 5x. A rejected step is retried with a smaller `h`.
- The solvers allocate a new state on every stage. That keeps the code short and is fine for small systems. A large system would want in-place stages.
- All methods are explicit, so stiff problems force `rk45` down to very small steps.
- Interpolation between points is linear. Dormand-Prince has a proper fourth-order dense output and this doesn't implement it, so take small steps if you need accurate values between points.

## 🚀 Quick start

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Compatible with Julia 1.9 and later. A quick taste, solving the oscillator above to `t = 10` and comparing with the closed form `(cos 10, -sin 10)`:

```
julia> last(sol)
2-element Vector{Float64}:
 -0.8390715287190943
  0.5440211106935285

julia> sol.accepted, sol.rejected
(234, 1)

julia> cos(10.0), -sin(10.0)
(-0.8390715290764524, 0.5440211108893698)
```

## 🧪 Testing

21 assertions across 10 test sets, all passing:

- Euler: accuracy against `exp(-2)` and the start and end points.
- `rk4`: accuracy, fourth-order convergence (halving the step cuts the error by a factor between 12 and 20), the final step landing on the end point, and a harmonic oscillator that returns to its start after a full period with energy held to 1e-8.
- `rk45`: meets its tolerance on the oscillator, adapts the step (rejections happen, step sizes vary, a looser tolerance takes fewer steps than a tighter one) and conserves the Lotka-Volterra invariant.
- Interpolation: a value between points, the end point and a `DomainError` outside the interval.
- Argument checks: a backwards `tspan` and a negative `dt` both throw `ArgumentError`.

## Licence

MIT.

