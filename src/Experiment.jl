using RigorousInvariantMeasures, LinearAlgebra, BallArithmetic
export Experiment, MultipleExperiments

function power_norms(A, N)
    norms = zeros(N)
    
    Aiter = A
    for i in 1:N
        norms[i] = BallArithmetic.svd_bound_L2_opnorm(Aiter)
        Aiter *= A
    end
    return norms
end

Γ(σ, K) = ((1 / (sqrt(σ^2 * 2 * π)))exp((-σ^2 * K^2 * π^2) / 2))

bound_ρ_σ_2(σ) = sqrt(1 / (sqrt(σ^2 * 2 * π)))

function deterministic_discretized(α, β, K)
    FFTNx = 8*K
    B = FourierAdjoint(K, FFTNx) 
    D(x) = T(x; α, β)
    PK = assemble(B, D)
    bPK = convert_matrix(PK)
    return bPK
end

function Experiment(α, β, σ, K; 
                        max_iter = 10, 
                        bPK = deterministic_discretized(α, β, K))
    
    bD = NoiseBall(σ, K)

    PσK = bD * bPK
    F = eigen(PσK.c)
    fσK = F.vectors[:, end] 
    fσK /= fσK[1]
    fσKs = symmetrize_density(fσK)
    bfσKs = BallVector(fσKs)
    residual = PσK * bfσKs - bfσKs
    ϵ = norm(residual.c, 2) + norm(residual.r, 2)
    @info "ϵ" ϵ
    lnn = FourierLogDer(α, β, K)

    λ = dot(lnn, fσKs)

    A = PσK[2:end, 2:end]
    norms = power_norms(A, max_iter)
    @info "Norms" norms

    valΓ = Γ(σ, K)
    @info "Γ" valΓ
    valρ = bound_ρ_σ_2(σ)
    @info "ρ" valρ

    coeff_err = ((1+valΓ+valρ)*valΓ+ϵ)
    err_L2 = (sum(norms)*coeff_err)/(1-norms[end])
    @info err_L2
    valΥ = Υ(α, β)
    @info valΥ
    return real(λ)+valΥ * hull(-err_L2, err_L2)
end

function MultipleExperiments(α, β, K, σ_arr)
    bPK = deterministic_discretized(α, β, K)
    lyap = [Experiment(α, β, interval(σ), K; bPK = bPK) for σ in σ_arr]
    return lyap
end