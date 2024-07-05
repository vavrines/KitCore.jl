sL = 0.13
sR = -0.27

KC.linear(sL, sR)
KC.vanleer(sL, sR)
KC.vanleer(sL, sR, sR)
KC.vanleer(sL, sR, sR, 3)
KC.minmod(sL, sR)
KC.minmod(sL, sR, sR)
KC.minmod(sL, sR, sR, 3)
KC.vanalbaba(sL, sR)
KC.superbee(0.5, 0.9)
KC.superbee(-1.0, -0.8)
KC.superbee(0.51, 2.01)
KC.weno5(-2.0, -1.0, 0.0, 1.0, 2.0)

#--- filter ---#
let deg = 2, u = rand(deg + 1)
    ℓ = rand(deg + 1)

    modal_filter!(u, 1e-6; filter = :l2)
    modal_filter!(u, 1e-6; filter = :l2opt)
    modal_filter!(u, 1e-6, ℓ; filter = :l1)
    modal_filter!(u, ℓ; filter = :lasso)
    modal_filter!(u, 10; filter = :exp)
    modal_filter!(u, 10; filter = :houli)
end

let deg = 2, u = rand(deg + 1, deg + 1)
    ℓ = rand(deg + 1, deg + 1)

    modal_filter!(u, 1e-6, 1e-6; filter = :l2)
    modal_filter!(u, 1e-6, 1e-6; filter = :l2opt)
    modal_filter!(u, 1e-6, 1e-6, ℓ; filter = :l1)
    modal_filter!(u, ℓ; filter = :lasso)
    modal_filter!(u, 2; filter = :exp)
end
