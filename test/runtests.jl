using Marl
using Test

decay(t, u) = -2u
oscillator(t, u) = [u[2], -u[1]]

@testset "Marl" begin
    @testset "euler" begin
        sol = euler(decay, 1.0, (0.0, 1.0); dt = 1e-4)
        @test isapprox(last(sol), exp(-2); atol = 1e-3)
        @test sol.t[1] == 0.0
        @test sol.t[end] ≈ 1.0
    end

    @testset "rk4 accuracy" begin
        sol = rk4(decay, 1.0, (0.0, 1.0); dt = 0.01)
        @test isapprox(last(sol), exp(-2); atol = 1e-8)
    end

    @testset "rk4 is fourth order" begin
        coarse = abs(last(rk4(decay, 1.0, (0.0, 1.0); dt = 0.1)) - exp(-2))
        fine = abs(last(rk4(decay, 1.0, (0.0, 1.0); dt = 0.05)) - exp(-2))
        @test 12 < coarse / fine < 20
    end

    @testset "final step lands on the end point" begin
        sol = rk4(decay, 1.0, (0.0, 1.0); dt = 0.3)
        @test sol.t[end] ≈ 1.0
        @test length(sol) == 5
    end

    @testset "systems" begin
        sol = rk4(oscillator, [1.0, 0.0], (0.0, 2pi); dt = 0.01)
        @test isapprox(last(sol), [1.0, 0.0]; atol = 1e-6)
        energy = [u[1]^2 + u[2]^2 for u in sol.u]
        @test maximum(abs.(energy .- 1)) < 1e-8
    end

    @testset "rk45 meets its tolerance" begin
        sol = rk45(oscillator, [1.0, 0.0], (0.0, 10.0); abstol = 1e-10, reltol = 1e-10)
        @test isapprox(last(sol), [cos(10.0), -sin(10.0)]; atol = 1e-7)
        @test sol.accepted == length(sol) - 1
    end

    @testset "rk45 adapts the step" begin
        stiffish(t, u) = -50 * (u - cos(t))
        sol = rk45(stiffish, 0.0, (0.0, 2.0); dt = 1.0)
        steps = diff(sol.t)
        @test sol.rejected > 0
        @test maximum(steps) > 2 * minimum(steps)
        loose = rk45(decay, 1.0, (0.0, 5.0); reltol = 1e-3, abstol = 1e-3)
        tight = rk45(decay, 1.0, (0.0, 5.0); reltol = 1e-9, abstol = 1e-9)
        @test loose.accepted < tight.accepted
        @test isapprox(last(tight), exp(-10); atol = 1e-8)
    end

    @testset "lotka-volterra stays on its invariant" begin
        lv(t, u) = [1.5u[1] - u[1] * u[2], -3u[2] + u[1] * u[2]]
        invariant(u) = u[1] - 3log(u[1]) + u[2] - 1.5log(u[2])
        sol = rk45(lv, [1.0, 1.0], (0.0, 20.0); abstol = 1e-10, reltol = 1e-10)
        @test isapprox(invariant(last(sol)), invariant([1.0, 1.0]); atol = 1e-6)
    end

    @testset "interpolation" begin
        sol = rk4(decay, 1.0, (0.0, 1.0); dt = 0.001)
        @test isapprox(sol(0.5005), exp(-1.001); atol = 1e-5)
        @test sol(1.0) == last(sol)
        @test_throws DomainError sol(2.0)
    end

    @testset "argument checks" begin
        @test_throws ArgumentError rk4(decay, 1.0, (1.0, 0.0); dt = 0.1)
        @test_throws ArgumentError euler(decay, 1.0, (0.0, 1.0); dt = -0.1)
    end
end

