using IntervalArithmetic

export aliasing_constant, aliasing_bound

@doc raw"""
    aliasing_constant(j, α, β) -> Interval

A constant ``C_j`` with ``|\hat f_j(n)| \le C_j/n^2`` for ``n \ne 0``, where
``f_j := e^{2\pi i j T}`` and ``T = T(\cdot; α, β)`` is the plateau map on ``[0,1]``
(`Dynamic.jl`), ``\hat f(n) := \int_0^1 f(y)e^{-2\pi i n y}\,dy``. Requires ``α > 1``.

Since ``T(0) = T(1) = 0``, ``f_j(0) = f_j(1) = 1``, and integrating by parts twice,
``|\hat f_j(n)| \le (|f_j'(1) - f_j'(0)| + \|f_j''\|_{L^1})/(4\pi^2 n^2)``, as in equation
`Pkell_decay` of Galatolo, Lopez Vereau, Marangio, Nisoli (the version on ``[-1,1]``).
With ``T'(y) = T_{-1,1}'(x)`` and ``T''(y) = 2T_{-1,1}''(x)``, ``x = 2y - 1``,
``|T_{-1,1}'(\pm1)| = α(1+β)``, ``\int_{-1}^1|T_{-1,1}''| = 2α(1+β)`` and
``\int_{-1}^1 T_{-1,1}'^2 = 2α^2(1+β)^2/(2α - 1)``, this gives
```math
C_j = \frac{2|j|α(1+β)}{\pi} + \frac{j^2α^2(1+β)^2}{2α - 1}.
```
"""
function aliasing_constant(j::Integer, α, β)
    inf(interval(α)) > 1 || throw(ArgumentError("need α > 1"))
    a = interval(α)
    b = interval(β)
    return 2 * abs(j) * a * (1 + b) / interval(π) + j^2 * a^2 * (1 + b)^2 / (2a - 1)
end

@doc raw"""
    aliasing_bound(j, K, N, α, β) -> Interval

Bound on ``|\hat f_j(\ell) - \hat f_{j,N}(\ell)|`` for ``|\ell| \le K``, where
``\hat f_{j,N}(\ell) := N^{-1}\sum_{m=0}^{N-1} f_j(m/N)e^{-2\pi i \ell m/N}`` is the discrete
coefficient from ``N`` samples and ``f_j = e^{2\pi i j T}``. Requires ``N \ge 4K``.

By the aliasing identity ``\hat f_{j,N}(\ell) = \sum_{q\in\mathbb Z}\hat f_j(\ell + qN)`` and
``|\hat f_j(n)| \le C_j/n^2`` (`aliasing_constant`), with ``|\ell + qN| \ge N(|q| - r)``,
``r := K/N``,
```math
\sum_{q \ne 0}|\hat f_j(\ell + qN)| \le \frac{2C_j}{N^2}\sum_{q\ge1}\frac{1}{(q - r)^2}
\le \frac{2C_j}{N^2}\Bigl(\frac{1}{(1 - r)^2} + \frac{1}{1 - r}\Bigr).
```
This replaces the constant ``\pi^2C_j/(3N^2)`` of Lemma `lem:Pkell_total_error` of the
paper above, whose proof uses ``|\ell + qN| \ge |q|N``, which fails for ``\ell \ne 0``.
"""
function aliasing_bound(j::Integer, K::Integer, N::Integer, α, β)
    N >= 4K || throw(ArgumentError("need N >= 4K"))
    r = interval(K) / interval(N)
    return 2 * aliasing_constant(j, α, β) / interval(N)^2 * (1 / (1 - r)^2 + 1 / (1 - r))
end
