
@everywhere function matrix_filename(α, β, K)
    αstr = replace(string(round(α, digits = 4)), "." => "p")
    βstr = replace(string(round(β, digits = 4)), "." => "p")
    return "unimodal_matrix_$(αstr)_$(βstr)_$(K).jld2"
end

@everywhere const local_matrix_cache = Dict{
    Tuple{Float64, Float64, Int}, Tuple{Float64, Any}}()

@everywhere function cleanup_matrix_cache!(max_age_seconds::Float64 = 600.0)
    now = time()
    keys_to_delete = [k
                      for (k, (ts, _)) in local_matrix_cache if now - ts > max_age_seconds]
    for k in keys_to_delete
        delete!(local_matrix_cache, k)
    end
    @debug "🧹 Cache cleanup: removed $(length(keys_to_delete)) matrices"
end

@everywhere using Logging


@everywhere function dowork(jobs, results)
    iteration = 0

    while true
        iteration += 1
        α, β, σ, K = take!(jobs)
        key = (α, β, K)

        λ = nothing
        norms = nothing
        t = @elapsed begin    
            with_logger(SimpleLogger(devnull, Logging.Error)) do
                # Suppress all worker output
                # Load or compute and update timestamp
                bPK = let
                    entry = get!(local_matrix_cache, key) do
                        (time(), PlateauExperiment.deterministic_discretized(α, β, K))
                    end
                    # Update access time
                    local_matrix_cache[key] = (time(), entry[2])
                    entry[2]
                end
                λ, norms = Experiment(α, β, σ, K; bPK = bPK)
            end
        end

        put!(results,
            (
                alpha = α,
                beta = β,
                sigma = σ,
                K = K,
                lambda = λ,
                norms = sup.(norms),
                time = t,
                id = myid()
            ))

        # Periodic cache cleanup every 100 iterations
        if iteration % 100 == 0
            cleanup_matrix_cache!()
        end
    end
end

function convergence_ok(res)::Bool
    λ = res.lambda
    return !(0 ∈ λ) || diam(λ) < 1e-13
end

function save_snapshot(basename, df, remaining, current_K, counter)
    suffix = (div(counter, 100) % 2 == 1) ? "a" : "b"
    fname = "$(basename)_snapshot_$(suffix).jld2"
    @debug "Saving snapshot to $fname (counter = $counter, remaining = $(length(remaining)))"
    JLD2.@save fname df remaining current_K counter
end

function resume_snapshot(basename, param_list, K0)
    df = nothing
    remaining = Set(param_list)
    fallback_current_K = Dict((α, β, σ) => K0 for (α, β, σ) in param_list)
    counter = 0

    # Check available snapshot files and sort by modification time
    files = filter(isfile, ["$(basename)_snapshot_a.jld2", "$(basename)_snapshot_b.jld2"])
    sorted_files = sort(files; by = x -> stat(x).mtime, rev = true)

    for file in sorted_files
        try
            JLD2.@load file df current_K counter

            # Compute which parameters are still missing
            computed_params = Set((row.alpha, row.beta, row.sigma) for row in eachrow(df))
            remaining = Set(p for p in param_list if p ∉ computed_params)

            # Ensure current_K contains all needed keys
            for (α, β, σ) in remaining
                if !haskey(current_K, (α, β, σ))
                    current_K[(α, β, σ)] = K0
                end
            end

            @info "✅ Loaded snapshot $(basename): $(length(remaining)) remaining."
            return df, remaining, current_K, counter
        catch e
            @warn "⚠️ Failed to load $file"
        end
    end

    # Fallback if no snapshot could be loaded
    @warn "⚠️ No valid snapshot could be loaded. Starting fresh."
    df = DataFrame(alpha = Float64[], beta = Float64[], sigma = Float64[],
        K = Int[], lambda = Interval[], norms = Vector{Float64}[], time = Float64[], id = Int[])
    return df, remaining, fallback_current_K, counter
end

function adaptive_dispatch_parallel(param_list, K0::Int,
        jobs, results; max_K = 512, basename = "results")
    df, remaining, current_K, counter = resume_snapshot(basename, param_list, K0)

    @debug "Remaining", length(remaining)

    # Submit one job per unresolved parameter
    for (α, β, σ) in remaining
        K = current_K[(α, β, σ)]
        @async put!(jobs, (α, β, σ, K))
    end

    @debug "Submitted jobs"
    while !isempty(remaining)
        @debug "Waiting on result... $(length(remaining)) remaining"
        res = take!(results)
        p = (res.alpha, res.beta, res.sigma)
        K = res.K

        if convergence_ok(res)
            if diam(res.lambda) < 1e-6
                push!(df, res)
                delete!(remaining, p)
                delete!(current_K, (p[1], p[2], p[3]))
                @debug "✓ Converged: $p at K=$K"
            else
                K′ = 2K
                if K′ > max_K
                    @warn "✗ Max resolution exceeded for $p"
                    push!(df, res)
                    delete!(remaining, p)
                    @info "Remaining $(length(remaining))"
                    delete!(current_K, (p[1], p[2], p[3]))
                else
                    current_K[(p[1], p[2], p[3])] = K′
                    #            ensure_matrix(p[1], p[2], K′)
                    @async put!(jobs, (p[1], p[2], p[3], K′))
                    @info "↻ λ too wide for $p. Retrying with K=$K′"
                end
            end
        else
            K′ = 2K
            if K′ > max_K
                @warn "✗ Gave up on $p due to 0 ∈ λ and K > $max_K"
                push!(df, res)
                delete!(remaining, p)
                @info "Remaining $(length(remaining))"
                delete!(current_K, (p[1], p[2], p[3]))
            else
                current_K[(p[1], p[2], p[3])] = K′
                @async put!(jobs, (p[1], p[2], p[3], K′))
                @debug "↻ 0 ∈ λ for $p. Retrying with K=$K′"
            end
        end

        counter += 1

        snap_freq = length(remaining) < 100 ? 10 : 100

        if counter % snap_freq == 0
            save_snapshot(basename, df, remaining, current_K, counter)
            @info "Progress: counter = $counter, remaining = $(length(remaining)), avg time/job = $(round(sum(df.time) / max(1, size(df, 1)), digits=4)) sec, est. time left = $(round(sum(df.time) / max(1, size(df, 1)) * length(remaining) / nworkers(), digits=2)) sec"
        end
    end

    save_snapshot(basename, df, remaining, current_K, counter)
    return df
end

function expand_norms(df::DataFrame; prefix = "norm_")
    L = length(df.norms[1])  # assume all have same length
    # Build new columns for each component
    for i in 1:L
        df[!, Symbol(prefix * string(i))] = getindex.(df.norms, i)
    end
    # Remove the original norms column
    select!(df, Not(:norms))
    return df
end