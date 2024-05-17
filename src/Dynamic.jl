
export T_minus_one_one, τ_1, τ_2, T, Υ

"""
Plateau map on [-1, 1]
"""
function T_minus_one_one(x; α, β)
    return β - (1 + β) * abs(x)^α
end

"""
Coordinate change from [-1,1] to [0,1]
"""
function τ_1(x)       # de [-1,1] a [0,1]
    return (x + 1) / 2
end

"""
Coordinate change from [0,1] to [-1,1]
"""
function τ_2(x)       # de [0,1] a [-1,1]
    return 2 * x - 1
end

"""
Plateau map on [0, 1]
"""
T(x; α, β) = τ_1(T_minus_one_one(τ_2(x); α = α, β = β))

@doc raw"""
L2-norm of log(|T'|).

This can be computed explictly as 

``\Upsilon(\alpha, \beta) = \sqrt{2}((\ln((1+\beta)\alpha)-(\alpha-1))^2+(\alpha-1)^2)^{\frac{1}{2}}``
"""
Υ(α, β) = sqrt(2) * ((log((β + 1) * α) - (α - 1))^2 + (α - 1)^2)^(1 / 2)
