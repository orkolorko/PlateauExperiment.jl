@everywhere function matrix_filename(α, β, K)
    αstr = replace(string(round(α, digits=4)), "." => "p")
    βstr = replace(string(round(β, digits=4)), "." => "p")
    return "unimodal_matrix_$(αstr)_$(βstr)_$(K).jld2"
end

@everywhere const local_matrix_cache = Dict{Tuple{Float64, Float64, Int}, Any}()

@everywhere function dowork(jobs, results)
    while true
        α, β, σ, K = take!(jobs)
        key = (α, β, K)
        bPK = get!(local_matrix_cache, key) do
            JLD2.load(matrix_filename(α, β, K))["P"]
        end
        t = @elapsed λ = Experiment(α, β, σ, K; bPK=bPK)
        
        put!(results, (
            alpha = α,
            beta = β,
            sigma = σ,
            K = K,
            lambda = λ,
            time = t,
            id = myid()
        ))
    end
end

function ensure_matrix(α, β, K::Int)
    fname = matrix_filename(α, β, K)
    if !isfile(fname)
        old_logger = current_logger()
        try
            # Suppress info/warn messages inside this block
            disable_logger = SimpleLogger(stderr, Logging.Error)
            with_logger(disable_logger) do
                P = PlateauExperiment.deterministic_discretized(α, β, K)
                JLD2.@save fname P
            end
        finally
            # Restore previous logger
            global_logger(old_logger)
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
    if isfile("$(basename)_snapshot_b.jld2")
        try
            JLD2.@load "$(basename)_snapshot_b.jld2" df remaining current_K counter
            return df, remaining, current_K, counter
        catch e
            @warn "Failed to load snapshot_b.jld2: $(e.message)"
        end
    end

    if isfile("$(basename)_snapshot_a.jld2")
        try
            JLD2.@load "$(basename)_snapshot_a.jld2" df remaining current_K counter
            return df, remaining, current_K, counter
        catch e
            @warn "Failed to load snapshot_a.jld2: $(e.message)"
        end
    end

    @warn "No valid snapshot could be loaded. Starting fresh."
    df = DataFrame(alpha=Float64[], beta=Float64[], sigma=Float64[],
                   K=Int[], lambda=Interval[], time=Float64[], id=Int[])
    remaining = Set(param_list)
    current_K = Dict(p => K0 for p in param_list)
    counter = 0
    return df, remaining, current_K, counter
end

function adaptive_dispatch_parallel(param_list, K0::Int,
                                    jobs, results; max_K=1024, basename="results")
    df, remaining, current_K, counter = resume_snapshot(basename, param_list, K0)
    
    @debug "Remaining", length(remaining)

    # Submit one job per unresolved parameter
    for (α, β, σ) in remaining
        K = current_K[(α, β, σ)]
        ensure_matrix(α, β, K)
        @async put!(jobs, (α, β, σ, K))
    end

    @debug "Submitted jobs"
    while !isempty(remaining)
        @debug "Waiting on result... $(length(remaining)) remaining"
        res = take!(results)
        p = (res.alpha, res.beta, res.sigma)
        K = res.K

        if convergence_ok(res)
            if diam(res.lambda) < 1e-10
                push!(df, res)
                delete!(remaining, p)
                @debug "✓ Converged: $p at K=$K"
            else
                K′ = 2K
                if K′ > max_K
                    @warn "✗ Max resolution exceeded for $p"
                    delete!(remaining, p)
                else
                    current_K[p] = K′
                    ensure_matrix(p[1], p[2], K′)
                    @async put!(jobs, (p[1], p[2], p[3], K′))
                    @debug "↻ λ too wide for $p. Retrying with K=$K′"
                end
            end
        else
            K′ = 2K
            if K′ > max_K
                @warn "✗ Gave up on $p due to 0 ∈ λ and K > $max_K"
                delete!(remaining, p)
            else
                current_K[p] = K′
                ensure_matrix(p[1], p[2], K′)
                @async put!(jobs, (p[1], p[2], p[3], K′))
                @debug "↻ 0 ∈ λ for $p. Retrying with K=$K′"
            end
        end

        counter += 1
        if counter % 100 == 0
            save_snapshot(basename, df, remaining, current_K, counter)
            @info "Progress: counter = $counter, remaining = $(length(remaining)), avg time/job = $(round(sum(df.time) / max(1, size(df, 1)), digits=4)) sec, est. time left = $(round(sum(df.time) / max(1, size(df, 1)) * length(remaining) / nworkers(), digits=2)) sec"
        end
    end

    save_snapshot(basename, df, remaining, current_K, counter)
    return df
end
