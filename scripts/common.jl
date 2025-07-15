
@everywhere function matrix_filename(α, β, K)
    return "unimodal_matrix_$(replace(string(α), '.' => 'p'))_$(replace(string(β), '.' => 'p'))_$(K).jld2"
end

@everywhere function dowork(jobs, results)
    while true
        job = take!(jobs)
        α, β, σ, K = job  # destructure the 4-tuple
        bPK = JLD2.load(matrix_filename(α, β, K))["P"]
        t = @elapsed λ = Experiment(α, β, σ, K; bPK = bPK)
        put!(results,
            (
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
        @info "Matrix for (α=$α, β=$β, K=$K) not found. Computing..."
        P = PlateauExperiment.deterministic_discretized(α, β, K)  # use dummy α, β if needed
        JLD2.@save fname P
    end
end

function convergence_ok(res)::Bool
    λ = res["lambda"]
    return !(0 ∈ λ) || diam(λ) < 1e-13
end

function adaptive_dispatch(param_list, K0, jobs, results)
    df = DataFrame(alpha = Float64[], beta = Float64[], sigma = Float64[],
        K = Int[], lambda = Interval[], time = Float64[], id = Int[])

    for (α, β, σ) in param_list
        K = K0
        while true
            ensure_matrix(α, β, K)
            put!(jobs, (α, β, σ, K))
            res = take!(results)
            if convergence_ok(res)
                if diam(res.lambda) < 1e-10
                    push!(df, res)
                    break
                else
                    @info "λ interval too wide (diam = $(diam(res.lambda))) for α=$α, σ=$σ, retrying with K=$(2K)"
                    K *= 2
                end
            end
        end
    end
    return df
end