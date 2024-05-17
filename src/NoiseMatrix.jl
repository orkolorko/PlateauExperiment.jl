using IntervalArithmetic, BallArithmetic

"""
Build a Matrix of Intervals representig Convolution with Gaussian Noise of variance σ
Galerkin truncated at frequence K
"""
function NoiseInterval(σ, K)
    return Diagonal([[exp((-σ^2 * π^2 * interval(k)^2) / 2) for k in 0:K];
                     [exp((-σ^2 * π^2 * interval(k)^2) / 2) for k in (-K):-1]])
end


"""
Build a Matrix of Balls representig Convolution with Gaussian Noise of variance σ
Galerkin truncated at frequence K
"""
function NoiseBall(σ, K)
    D = NoiseInterval(σ, K)
    Dc, Dr = IntervalArithmetic.mid.(D), IntervalArithmetic.radius.(D)
    bD = BallMatrix(Dc, Dr)
end