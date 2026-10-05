using IntervalArithmetic, BallArithmetic

# Both matrices enclose the same Galerkin matrix, so they must intersect entrywise; with N = 4K
# the aliasing term is the dominant part of the radius, so this fails if aliasing_bound is too small.
let α = 3.5, β = 1.0, K = 16
    P1 = deterministic_discretized(α, β, K; N = 4K)
    P2 = deterministic_discretized(α, β, K; N = 2^16)
    d = abs.(P1.c .- P2.c)
    @test all(d .<= P1.r .+ P2.r)
    @test maximum(P1.r) > 100 * maximum(P2.r)
end

@test_throws ArgumentError aliasing_bound(1, 16, 32, 3.5, 1.0)
@test_throws ArgumentError aliasing_constant(1, 1.0, 1.0)
