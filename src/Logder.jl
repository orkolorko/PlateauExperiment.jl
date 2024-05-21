f(t) = sin(t) / t

function coeff_eval_at_0(i, x)
    return x^(2i + 1) / (factorial(big(2i + 1)) * (2i + 1))
end

using TaylorModels

function integrate_i_i_plus_1(f, i; degree)            #aqui puede usar porque evalua las derivadas em m que es diferente de 0
    I = (@interval i*π (i + 1)*π)                             #esto es el intervalo [i*π,(i+1)*π], i>0
    m = Interval(IntervalArithmetic.mid(I))                                    #mid(I) ponto medio do intervalo I, Interval(mid(I))=un intervalo pequenho que contenga mid(I)
    tm = TaylorModel1(degree, m, I)                         #ese m es el valor que será evaluada en la serie de Taylos i.e f(m)+f'(m)t/1!+...+f^{(400)}(m)t^{400}/400!+erro
    prim = TaylorSeries.integrate(f(tm))                    #esto es la primitiva de f, al poner f(tm), recien esta diciendo quien va ser la función f para la cual voy hacer el polinomi de taylor i.e. para obtener la sumatoria que esta en la linea de arriba
    return prim(I.hi - m) - prim(I.lo - m)                        #I.hi, I.lo parte superior e inferior do intervalo I respectivamente, como m é um intervalo, vai ficar o intevalo prim([I.hi-m1,I.hi-m2])-prim([I.lo-m1,I.lo-m2]) where denotamos m=[m1,m2]
end

@doc raw"""
Compute a vector of K elements containing
```math
\frac{1}{j\pi}\int_0^{j\pi} \frac{sin(t)}{t}dt
```
for ``1 \leq j \leq K`` by using Taylor Models integration for ``[\pi, j\pi]``.
"""
function IntegralSintOvert(K; N = 1000, degree = 40)
    Pi = @biginterval π

    # We start by computing the integral of sin(t)/t on [0, π]                              
    v = [(-1)^i * coeff_eval_at_0(i, Pi) for i in 0:N]
    # Since it is an alternating sum, the error in the sum is
    error = Interval(-1, 1) * abs(coeff_eval_at_0(N + 1, Pi))
    # the following is an enclosure of the integral of sin(t)/t on [0, π]
    integral_with_error = sum(v) + error
    val = [integral_with_error;
           [integrate_i_i_plus_1(f, i; degree = degree) for i in 1:(K - 1)]]
    cum_sum = cumsum(val)

    #The vector `coeff` contains at the index $i$ the value of
    # \frac{1}{i\pi}\int_0^{i\pi} \frac{\sin(t)}{t}dt.
    coeff = cum_sum ./ ([i * @interval π for i in 1:K])

    return coeff
end

function FourierLogDer(α, β, K; N = 1000, degree = 40,
        coeff = IntegralSintOvert(K; N = N, degree = degree))
    coeff_0_1 = [(-1)^i for i in 1:length(coeff)] .* coeff
    coeff_fft = [coeff_0_1; reverse(coeff_0_1[1:end])]
    ln = [log((1 + β) * α) - (α - 1); -(α - 1) * coeff_fft]
    # we need to translate from Fourier coefficients on [-1,1] 
    # to Fourier coefficients on [0,1]
    return ln
end