import Pkg;
Pkg.activate(@__DIR__)
#Pkg.instantiate()

using Logging, Dates, Distributed, LinearAlgebra, SlurmClusterManager, DataFrames, JLD2, CSV

datetime = Dates.now()

global_t0 = time()

if haskey(ENV, "SLURM_NTASKS")
    procs = addprocs(SlurmManager())
    location = "slurm"
else
    procs = addprocs(4)
    location = "local"
end

nprocs = length(procs)

@everywhere import Pkg
@everywhere Pkg.activate(@__DIR__)
#@everywhere Pkg.instantiate()

@everywhere using LinearAlgebra, PlateauExperiment, JLD2, IntervalArithmetic

const job_channel = RemoteChannel(() -> Channel{Tuple{Float64, Float64, Float64, Int64}}(1024))
const result_channel = RemoteChannel(() -> Channel{NamedTuple}(1024))

include("./common.jl")

foreach(
    pid -> remote_do(dowork, pid, job_channel, result_channel),
    workers()
)

sigma_0 = 1/16

N = 512

param_list = [(3.0, β, σ) for β in range(start = 51/64, length = N, step = 1/(8*N)), σ in range(start = sigma_0, length = N, step = (1-sigma_0)/N)]

df = adaptive_dispatch_parallel(param_list, 64, job_channel, result_channel, basename = "NIC")

# Add derived columns
df.lambda_lo = inf.(df.lambda)
df.lambda_hi = sup.(df.lambda)

# Save CSV (omit original lambda if you prefer)
CSV.write("results_NIC.csv", select(df, Not(:lambda)))

# Save full object (intervals intact)
JLD2.@save "results_NIC.jld2" df

max_diam = maximum(diam.(df.lambda))
elapsed_time = time() - global_t0

@info "Maximum λ diameter: $max_diam"
@info "Total elapsed time: $(round(elapsed_time; digits=2)) seconds"


rmprocs(procs)