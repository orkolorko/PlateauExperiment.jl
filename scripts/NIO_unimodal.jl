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

sigma_0 = 1/128+1/512+1/8192

param_list = [(α, 1.0, σ) for α in range(start = 3.0, length = 1024, step = 1/1024), σ in range(start = sigma_0, length = 1024, step = (1-sigma_0)/1024)]

df = adaptive_dispatch_parallel(param_list, 64, job_channel, result_channel)

# Add derived columns
df.lambda_lo = inf.(df.lambda)
df.lambda_hi = sup.(df.lambda)

# Save CSV (omit original lambda if you prefer)
CSV.write("results_NIO.csv", select(df, Not(:lambda)))

# Save full object (intervals intact)
JLD2.@save "results_NIO.jld2" df

max_diam = maximum(diam.(df.lambda))
elapsed_time = time() - global_t0

@info "Maximum λ diameter: $max_diam"
@info "Total elapsed time: $(round(elapsed_time; digits=2)) seconds"


rmprocs(procs)