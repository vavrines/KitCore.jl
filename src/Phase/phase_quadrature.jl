# ============================================================
# Quadrature Methods
# ============================================================

"""
$(SIGNATURES)

Maxwell quadrature

## Arguments
* `N`: quadrature order (MUST less than 33)
"""
function maxwell_quadrature(N2, C = 1.0)
    N = N2 ÷ 2

    a = zeros(N)
    b = zeros(N)
    a[1] = 1.0 / sqrt(pi)
    a[2] = 2.0 / sqrt(pi) / (pi - 2.0)
    b[2] = a[1] / (a[1] + a[2]) / 2.0

    for i = 3:N
        b[i] = (i - 2) + 1.0 / 2.0 - b[i-1] - a[i-1]^2
        a[i] = ((i - 1)^2 / (4.0 * b[i]) - b[i-1] - 1.0 / 2) / a[i-1] - a[i-1]
    end

    J = Diagonal(a) + diagm(1 => sqrt.(b[2:N])) + diagm(-1 => sqrt.(b[2:N]))

    v, V = eigen(J)

    w = V[1, :] .^ 2 .* sqrt(pi) / 2.0

    vw = [v w]
    vw = sortslices(vw, dims = 1, by = x -> x[1])
    v = vw[:, 1]
    w = vw[:, 2]

    Xis = vcat(-reverse(v), v)
    weights = vcat(reverse(w), w)
    weights .*= exp.(Xis .^ 2) .* C
    Xis .*= C

    return Xis, weights
end
