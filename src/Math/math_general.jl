# ============================================================
# General Mathematical Methods
# ============================================================

"""
$(SIGNATURES)

Python linspace function
"""
linspace(start, stop, n::Integer) = collect(range(start, stop = stop, length = n))


"""
$(SIGNATURES)

Heaviside step function
"""
heaviside(x) = ifelse(x >= 0, one(x), zero(x))


"""
$(SIGNATURES)

Fortran sign function
"""
fortsign(x, y) = abs(x) * sign(y)


"""
$(SIGNATURES)

Calculate rate of convergence
"""
convergence_order(e1, e2, r = 2) = log(e1 / e2) / log(r)


"""
$(SIGNATURES)

Calculate L1 error
"""
L1_error(u::T, ue::T, Δx) where {T<:AA} = sum(abs.(u .- ue) .* Δx)


"""
$(SIGNATURES)

Calculate L2 error
"""
L2_error(u::T, ue::T, Δx) where {T<:AA} = sqrt(sum((abs.(u .- ue) .* Δx) .^ 2))


"""
$(SIGNATURES)

Calculate L∞ error
"""
L∞_error(u::T, ue::T, Δx) where {T<:AA} = maximum(abs.(u .- ue) .* Δx)


"""
$(SIGNATURES)

Evaluate quadrature weight from Newton-Cotes rule
"""
function newton_cotes(idx::T, num::T) where {T<:Integer}
    if idx == 1 || idx == num
        nc_coeff = 14.0 / 45.0
    elseif (idx - 5) % 4 == 0
        nc_coeff = 28.0 / 45.0
    elseif (idx - 3) % 4 == 0
        nc_coeff = 24.0 / 45.0
    else
        nc_coeff = 64.0 / 45.0
    end

    return nc_coeff
end


"""
$(SIGNATURES)

Gauss Legendre integral for fast spectral method

## Arguments
* `N`: number of quadrature points, 
* `a, b`: integral range
"""
function lgwt(N::Integer, a::Real, b::Real)
    x = zeros(N)
    w = zeros(N)

    N1 = N
    N2 = N + 1

    y = zeros(N1)
    y0 = zeros(N1)
    Lp = zeros(N1)
    L = zeros(N1, N2)

    # initial guess
    for i = 1:N1
        y[i] =
            cos((2.0 * (i - 1.0) + 1.0) * 4.0 * atan(1.0) / (2.0 * (N - 1.0) + 2.0)) +
            0.27 / N1 *
            sin(4.0 * atan(1.0) * (-1.0 + i * 2.0 / (N1 - 1.0)) * (N - 1.0) / N2)
        y0[i] = 2.0
    end

    # compute the zeros of the N+1 legendre Polynomial
    # using the recursion relation and the Newton method
    while maximum(abs.(y .- y0)) > 0.0000000000001
        L[:, 1] .= 1.0
        L[:, 2] .= y
        for k = 2:N1
            @. L[:, k+1] = ((2.0 * k - 1.0) * y * L[:, k] - (k - 1) * L[:, k-1]) / k
        end
        @. Lp = N2 * (L[:, N1] - y * L[:, N2]) / (1.0 - y^2)
        @. y0 = y
        @. y = y0 - L[:, N2] / Lp
    end

    # linear map from [-1 1] to [a,b]
    @. x = (a * (1.0 - y) + b * (1.0 + y)) / 2.0
    @. w = N2^2 * (b - a) / ((1.0 - y^2) * Lp^2) / N1^2

    return x, w
end
