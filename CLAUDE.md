# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**PlateauExperiment.jl** is a Julia package for rigorous numerical study of noise-induced chaos (NIC) and noise-induced oscillations (NIO) in plateau maps. It computes Lyapunov exponents with certified error bounds using interval/ball arithmetic, Fourier analysis, and distributed parallel processing.

Authors: Isaia Nisoli, Charles Lopes Vereau.

## Commands

### Run Tests
```bash
# From repo root
julia --project=. -e "using Pkg; Pkg.test()"

# Run a single test file interactively
julia --project=. test/TestExperiment.jl
```

### Install Dependencies
```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
# Scripts have their own environment:
julia --project=scripts/ -e "using Pkg; Pkg.instantiate()"
```

### Run Experiments Locally
```bash
cd scripts/
julia --project=. NIO_unimodal.jl 100   # 100×100 parameter grid
julia --project=. NIC_unimodal.jl 100
```

### Run on HPC Cluster (SLURM)
```bash
cd scripts/
sbatch NIO.sbatch   # 32 tasks, 72hr limit
sbatch NIC.sbatch
```

## Architecture

### Package Structure (`src/`)

- **`Dynamic.jl`** — Plateau map `T(x; α, β)` on [0,1] and coordinate transforms `τ_1`, `τ_2`. The map `T_minus_one_one(x; α, β) = β - (1+β)|x|^α` is the canonical form on [-1,1]. `Υ(α, β)` computes the L2-norm of `log|T'|`.

- **`Operators.jl`** — Builds Gaussian noise operators: `NoiseInterval(σ, K)` returns a diagonal matrix of interval Fourier coefficients; `NoiseBall(σ, K)` converts to ball arithmetic. Also handles `convert_matrix`, `symmetrize_density`, and `compute_residual`.

- **`Logder.jl`** — Rigorously computes Fourier coefficients of `log|T'|` via Taylor models (`TaylorModels.jl`). `FourierLogDer(α, β, K)` is the main entry point.

- **`Experiment.jl`** — Core computational pipeline:
  1. `deterministic_discretized(α, β, K)` — Fourier-discretized transfer operator using `RigorousInvariantMeasures`
  2. `power_norms(A, N)` — Computes `‖A^i‖` for i=1..N via `safe_svd_bound`
  3. `Experiment(α, β, σ, K)` — Full pipeline: builds `P_σ^K = B_σ * P^K`, finds invariant density eigenvector, computes Lyapunov exponent as `⟨density, log|T'|⟩`, returns `(λ::Interval, norms::Vector{Interval})`
  4. `MultipleExperiments(α, β, K, σ_arr)` — Batch version

### Scripts (`scripts/`)

- **`common.jl`** — Parallel job orchestration:
  - `dowork(jobs, results)` — Worker loop: pulls `(α, β, σ, K)` from channel, caches deterministic matrices by filename key, runs `Experiment`, pushes results
  - `convergence_ok(res)` — Checks validity of λ interval (non-missing, reasonable diameter)
  - `adaptive_dispatch_parallel(...)` — Main orchestrator: manages job/result channels, implements adaptive K doubling on convergence failure, checkpoints every 100 iterations to alternating `.jld2` snapshot files

- **`NIO_unimodal.jl`** — NIO parameter sweep: σ ∈ [1/128, 1/16], α ∈ [3.0, 4.0], β=1.0 fixed, K starts at 64
- **`NIC_unimodal.jl`** — NIC parameter sweep: σ ∈ [1/16, 1], β ∈ [51/64, 63/64], α=3.0 fixed, K starts at 64

### Key Design Patterns

- **Rigorous arithmetic:** All core computations use `IntervalArithmetic` / `BallArithmetic` to certify error bounds. Never replace interval types with floats.
- **Adaptive K refinement:** When `convergence_ok` fails, K is doubled and the job is re-queued. This is the primary mechanism for ensuring accuracy.
- **Matrix caching:** Workers cache the expensive `deterministic_discretized` computation by `(α, β, K)` key; `cleanup_matrix_cache!` runs periodic LRU eviction.
- **Fault-tolerant snapshots:** `adaptive_dispatch_parallel` alternates saving to `_a.jld2` / `_b.jld2` snapshots; prior results can be loaded to resume interrupted runs.
- **Two output formats:** JLD2 preserves interval types for further computation; CSV expands norm vectors to individual columns for analysis.

### Dependencies
- `RigorousInvariantMeasures` — Fourier basis operators and `FourierAdjoint`
- `IntervalArithmetic`, `BallArithmetic` — Certified numerics
- `TaylorModels` — Rigorous integration in `Logder.jl`
- `FFTW` — Fourier transforms
- `Distributed`, `SlurmClusterManager` — Parallelism
- `JLD2`, `CSV`, `DataFrames` — Output
