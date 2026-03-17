#!/usr/bin/env julia
# Script to recompute Table 1 (tab:valNIO), Table 2 (tab:valNIC), and Table 3 (tab:crossing)
# with updated code: 256*K oversampling + corrected Theorem 6.5 error bounds

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PlateauExperiment
using IntervalArithmetic

# ── Helper ────────────────────────────────────────────────────────────────────
function run_point(α, β, σ, K)
    bPK = PlateauExperiment.deterministic_discretized(α, β, K)
    λ, _ = Experiment(α, β, interval(σ), K; bPK = bPK, max_iter = 20)
    return λ
end

# ── Table 1: NIO (β = 1.0) ────────────────────────────────────────────────────
println("\n=== TABLE 1: tab:valNIO (β=1.0, K=128) ===")
flush(stdout)

NIO_points = [
    (3.25, 1.0, 0.2),
    (3.5,  1.0, 0.2),
    (3.5,  1.0, 0.4),
    (3.8,  1.0, 0.4),
]

for (α, β, σ) in NIO_points
    λ = run_point(α, β, σ, 128)
    println("α=$α, β=$β, σ=$σ  →  λ = $λ   (diam=$(diam(λ)))")
    flush(stdout)
end

# ── Table 2: NIC (α = 3.0) ────────────────────────────────────────────────────
println("\n=== TABLE 2: tab:valNIC (α=3.0, K=128) ===")
flush(stdout)

NIC_points = [
    (3.0, 0.875, 0.2),
    (3.0, 0.9,   0.2),
    (3.0, 0.875, 0.4),
    (3.0, 0.9,   0.4),
]

for (α, β, σ) in NIC_points
    λ = run_point(α, β, σ, 128)
    println("α=$α, β=$β, σ=$σ  →  λ = $λ   (diam=$(diam(λ)))")
    flush(stdout)
end

# ── Table 3: Crossings (K=256) ────────────────────────────────────────────────
println("\n=== TABLE 3: tab:crossing (K=256) ===")
flush(stdout)

crossing_points = [
    (1,  3.0,      0.523926, 0.525757),
    (2,  3.04883,  0.498291, 0.500122),
    (3,  3.09766,  0.474487, 0.476318),
    (4,  3.14648,  0.452515, 0.454346),
    (5,  3.19531,  0.430542, 0.432373),
    (6,  3.24414,  0.4104,   0.412231),
    (7,  3.29297,  0.390259, 0.39209),
    (8,  3.3418,   0.370117, 0.371948),
    (9,  3.39062,  0.349976, 0.351807),
    (10, 3.43945,  0.328003, 0.329834),
    (11, 3.48828,  0.302368, 0.304199),
    (12, 3.53711,  0.273071, 0.274902),
    (13, 3.58594,  0.23645,  0.238281),
    (14, 3.63477,  0.194336, 0.196167),
    (15, 3.68359,  0.165039, 0.16687),
    (16, 3.73242,  0.146729, 0.14856),
    (17, 3.78125,  0.133911, 0.135742),
    (18, 3.83008,  0.124756, 0.126587),
    (19, 3.87891,  0.115601, 0.117432),
    (20, 3.92773,  0.110107, 0.111938),
    (21, 3.97656,  0.104614, 0.106445),
]

for (row, α, σ1, σ2) in crossing_points
    β = 1.0
    bPK = PlateauExperiment.deterministic_discretized(α, β, 256)
    λ1, _ = Experiment(α, β, interval(σ1), 256; bPK = bPK, max_iter = 20)
    λ2, _ = Experiment(α, β, interval(σ2), 256; bPK = bPK, max_iter = 20)
    println("Row $row: α=$α")
    println("  σ1=$σ1: λ1 = $λ1")
    println("  σ2=$σ2: λ2 = $λ2")
    flush(stdout)
end

println("\n=== DONE ===")
