
"""
Builds a Matrix representig Convolution with Gaussian Noise of variance σ
Galerkin truncated at frequence K
"""
function D(σ, K)
    return Diagonal([[exp((-σ^2 * π^2 * interval(k)^2) / 2) for k in 0:K];
                     [exp((-σ^2 * π^2 * interval(k)^2) / 2) for k in (-K):-1]])
end