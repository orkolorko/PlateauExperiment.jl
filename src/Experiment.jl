using RigorousInvariantMeasures, LinearAlgebra, BallArithmetic
using FFTW   # loads the FFTWExt extension of RigorousInvariantMeasures, which assembles FourierAdjoint
export Experiment, MultipleExperiments, deterministic_discretized


function safe_svd_bound(A)
    # Step 1: Try fast upper bound
    ub = BallArithmetic.upper_bound_L2_opnorm(A)
    if ub < 1e-2
        @debug "✅ SVD skipped: upper bound is small (ub = $ub)"
        return ub  
    end

    # Step 2: Attempt full SVD, fallback if needed
    try
        return BallArithmetic.svd_bound_L2_opnorm(A)  # or your svd_bound_L2_opnorm(A)
    catch e
        @warn "⚠️ SVD failed, falling back to upper_bound_L2_opnorm: $e"
        return ub
    end
end

function power_norms(A, N)
    norms = zeros(N)
    
    K = N

    Aiter = A
    for i in 1:N

        norms[i] = safe_svd_bound(Aiter)
        Aiter *= A
        # if norms[i]<1
        #     K = i
        #     break
        # end
    end

    return norms
end

# L¹→L² tail bound Γ_{σ,K} (Theorem 6.3 / Theorem on tail estimate)
Γ(σ, K) = sqrt(coth(1 / (2 * σ^2)) / (σ * sqrt(π))) * exp((-σ^2 * K^2 * π^2) / 2)

# L¹→L¹ tail bound Γ^(1)_{σ,K} (Lemma 8.7, explicit integral bound)
Γ1(σ, K) = (2 / (σ^2 * π^2 * K)) * exp((-σ^2 * K^2 * π^2) / 2)

# ‖ρ_σ‖_{L²(ℝ)} · √coth(1/2σ²)  (enters second term of δ in Theorem 6.5)
bound_ρ_coth(σ) = sqrt(1 / (2 * σ * sqrt(π))) * sqrt(coth(1 / (2 * σ^2)))


# function process_norms(norms)
#     rest = 0
#     contract = false
#     for i in 1:length(norms)
#         if norms[i] < 1
#             contract = true
#         end
#         if norms[i] > 1 && contract == true
#             rest = i
#             break
#         end
#     end
#     good_norms = norms[1:rest]
#     return good_norms
# end


@doc raw"""
    deterministic_discretized(α, β, K; N = 1024K) -> BallMatrix

Enclosure of the Fourier--Galerkin matrix ``P_{jk} = \int_0^1 e^{2\pi i k y}e^{-2\pi i jT(y)}\,dy``,
``|j|, |k| \le K``, in the layout ``[0:K; -K:-1]``. RigorousInvariantMeasures computes the discrete
coefficients from ``N`` samples, with the map evaluated on interval sample points and the FFT
rounding enclosed; the aliasing error of row ``j`` is added here, by `aliasing_bound`.
"""
function deterministic_discretized(α, β, K; N = 1024 * K)
    B = FourierAdjoint(K, N)
    D(x) = T(x; α, β)
    PK = assemble(B, D)
    for r in axes(PK, 1)
        j = r <= K + 1 ? r - 1 : r - 1 - (2K + 1)
        A = aliasing_bound(j, K, N, α, β)
        e = interval(-sup(A), sup(A))
        PK[r, :] .= PK[r, :] .+ complex(e, e)
    end
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
    @debug "ϵ" ϵ
    lnn = FourierLogDer(α, β, K)

    λ = dot(lnn, fσKs)

    A = PσK[2:end, 2:end]
    norms = interval.(power_norms(A, max_iter))

    
    @debug "Norms" 
    @debug norms[1]
    @debug norms[end]

    valΓ  = Γ(σ, K)
    valΓ1 = Γ1(σ, K)
    valρc = bound_ρ_coth(σ)
    @debug "Γ"  valΓ
    @debug "Γ1" valΓ1
    @debug "ρ·√coth" valρc

    # δ from Theorem 6.5: Γ(1+Γ¹) + ‖ρ‖·√coth·Γ¹, with ‖f_σ‖_L¹ = 1
    δ = valΓ * (1 + valΓ1) + valρc * valΓ1
    coeff_err = (δ + ϵ)
    err_L2 = (sum(norms)*coeff_err)/(1-norms[end])
    @debug "err_L2" err_L2
    valΥ = Υ(α, β)
    @debug "Υ" valΥ
    @debug "diam λ" diam(real(λ))

    return real(λ)+valΥ * err_L2 * interval(-1,1), norms
end

function MultipleExperiments(α, β, K, σ_arr)
    bPK = deterministic_discretized(α, β, K)
    lyap = [Experiment(α, β, interval(σ), K; bPK = bPK) for σ in σ_arr]
    return lyap
end