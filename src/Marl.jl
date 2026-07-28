module Marl

export Solution, euler, rk4, rk45

struct Solution{T}
    t::Vector{Float64}
    u::Vector{T}
    accepted::Int
    rejected::Int
end

Base.length(sol::Solution) = length(sol.t)
Base.last(sol::Solution) = sol.u[end]

function (sol::Solution)(t::Real)
    sol.t[1] <= t <= sol.t[end] || throw(DomainError(t, "outside the solved interval"))
    i = searchsortedlast(sol.t, t)
    i == length(sol.t) && return sol.u[end]
    w = (t - sol.t[i]) / (sol.t[i + 1] - sol.t[i])
    return (1 - w) .* sol.u[i] .+ w .* sol.u[i + 1]
end

function span(tspan)
    t0, t1 = Float64(tspan[1]), Float64(tspan[2])
    t1 > t0 || throw(ArgumentError("tspan must run forwards"))
    return t0, t1
end

euler_step(f, t, u, h) = u .+ h .* f(t, u)

function rk4_step(f, t, u, h)
    k1 = f(t, u)
    k2 = f(t + h / 2, u .+ (h / 2) .* k1)
    k3 = f(t + h / 2, u .+ (h / 2) .* k2)
    k4 = f(t + h, u .+ h .* k3)
    return u .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

function fixed(step, f, u0, tspan, dt)
    t0, t1 = span(tspan)
    dt > 0 || throw(ArgumentError("dt must be positive"))
    t = t0
    u = float.(u0)
    ts = [t]
    us = [u]
    while t1 - t > 1e-12 * max(1.0, abs(t1))
        h = min(dt, t1 - t)
        u = step(f, t, u, h)
        t += h
        push!(ts, t)
        push!(us, u)
    end
    return Solution(ts, us, length(ts) - 1, 0)
end

euler(f, u0, tspan; dt) = fixed(euler_step, f, u0, tspan, dt)
rk4(f, u0, tspan; dt) = fixed(rk4_step, f, u0, tspan, dt)

const A = (
    (1 / 5,),
    (3 / 40, 9 / 40),
    (44 / 45, -56 / 15, 32 / 9),
    (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
    (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
    (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84),
)
const C = (1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0)
const B4 = (5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40)

function error_norm(e, u, unew, abstol, reltol)
    scaled = e ./ (abstol .+ reltol .* max.(abs.(u), abs.(unew)))
    return sqrt(sum(abs2, scaled) / length(scaled))
end

function rk45(f, u0, tspan; abstol = 1e-8, reltol = 1e-6, dt = nothing, maxsteps = 100_000)
    t0, t1 = span(tspan)
    t = t0
    u = float.(u0)
    h = dt === nothing ? (t1 - t0) / 100 : Float64(dt)
    ts = [t]
    us = [u]
    accepted = 0
    rejected = 0
    k1 = f(t, u)

    while t1 - t > 1e-12 * max(1.0, abs(t1))
        accepted + rejected < maxsteps || error("rk45 gave up after $maxsteps steps at t = $t")
        h = min(h, t1 - t)
        ks = Any[k1]
        for stage in 1:6
            increment = zero(u)
            for j in 1:stage
                increment = increment .+ A[stage][j] .* ks[j]
            end
            push!(ks, f(t + C[stage] * h, u .+ h .* increment))
        end

        fifth = zero(u)
        fourth = zero(u)
        for j in 1:6
            fifth = fifth .+ A[6][j] .* ks[j]
        end
        for j in 1:7
            fourth = fourth .+ B4[j] .* ks[j]
        end
        unew = u .+ h .* fifth
        err = error_norm(h .* (fifth .- fourth), u, unew, abstol, reltol)

        if err <= 1
            t += h
            u = unew
            k1 = ks[7]
            push!(ts, t)
            push!(us, u)
            accepted += 1
        else
            rejected += 1
        end

        factor = err == 0 ? 5.0 : clamp(0.9 * err^(-1 / 5), 0.2, 5.0)
        h *= factor
    end

    return Solution(ts, us, accepted, rejected)
end

end
