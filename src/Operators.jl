using IntervalArithmetic, BallArithmetic, LinearAlgebra

export NoiseInterval, NoiseBall, convert_matrix, symmetrize_density, compute_residual

"""
Build a Matrix of Intervals representig Convolution with Gaussian Noise of variance σ
Galerkin truncated at frequence K
"""
function NoiseInterval(σ, K)
    return Diagonal([[exp((-σ^2 * π^2 * interval(k)^2) / 2) for k in 0:K];
                     [exp((-σ^2 * π^2 * interval(k)^2) / 2) for k in (-K):-1]])
end

"""
Convert a matrix of Complex{Interval} to a BallMatrix with a center matrix of complex numbers
"""
function convert_matrix(M)
    center = IntervalArithmetic.mid.(real.(M)) + im * IntervalArithmetic.mid.(imag.(M))
    radius = sqrt.(IntervalArithmetic.radius.(real.(M)) .^ 2 +
                   IntervalArithmetic.radius.(imag.(M)) .^ 2)
    bM = BallMatrix(center, radius)
    return bM
end

"""
Build a Matrix of Balls representig Convolution with Gaussian Noise of variance σ
Galerkin truncated at frequence K
"""
NoiseBall(σ, K) = convert_matrix(NoiseInterval(σ, K))

"""
Symmetrize density
"""
function symmetrize_density(v)
    w = zeros(eltype(v), length(v))
    N = (length(v) - 1) ÷ 2

    w[1:(N + 1)] = v[1:(N + 1)]
    w[(end - N + 1):end] = [x' for x in reverse(v[2:(N + 1)])]

    return w
end

"""
Compute 
`Pv-v`
"""
function compute_residual(P, v)
    return P*v-v
end