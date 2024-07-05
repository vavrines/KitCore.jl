#--- 1D ---#
prim = zeros(3, 2)
prim[:, 1] .= [1.0, 3.0, 2.0]
prim[:, 2] .= [1 / 1864, -0.3, 1 / 1864]
# pure
Mu1, Mxi1, MuL1, MuR1 = KitCore.gauss_moments(prim[:, 1], 2)
Mu2, Mxi2, MuL2, MuR2 = KitCore.gauss_moments(prim[:, 2], 2)
# mixture
Mu, Mxi, MuL, MuR = KitCore.mixture_gauss_moments(prim, 2)

@test Mu[:, 1] ≈ Mu1 atol = 0.01
@test Mxi[:, 1] ≈ Mxi1 atol = 0.01
@test MuL[:, 1] ≈ MuL1 atol = 0.01
@test MuR[:, 1] ≈ MuR1 atol = 0.01
@test Mu[:, 2] ≈ Mu2 atol = 0.01
@test Mxi[:, 2] ≈ Mxi2 atol = 0.01
@test MuL[:, 2] ≈ MuL2 atol = 0.01
@test MuR[:, 2] ≈ MuR2 atol = 0.01

#--- 2D ---#
prim = zeros(4, 2)
prim[:, 1] .= [1.0, 3.0, 0.5, 1.3]
prim[:, 2] .= [1 / 1864, -1.3, 0.88, 2 / 1864]
# pure
KitCore.gauss_moments(prim[:, 1])
Mu1, Mv1, Mxi1, MuL1, MuR1 = KitCore.gauss_moments(prim[:, 1], 1)
Mu2, Mv2, Mxi2, MuL2, MuR2 = KitCore.gauss_moments(prim[:, 2], 1)
# mixture
Mu, Mv, Mxi, MuL, MuR = KitCore.mixture_gauss_moments(prim, 1)

Muv1 = KitCore.moments_conserve(MuR1, Mv1, Mxi1, 3, 1, 1)
Muv2 = KitCore.moments_conserve(MuR2, Mv2, Mxi2, 3, 1, 1)
Muv = KitCore.mixture_moments_conserve(MuR, Mv, Mxi, 3, 1, 1)

@test Muv[:, 1] ≈ Muv1 atol = 0.01
@test Muv[:, 2] ≈ Muv2 atol = 0.01

#--- 3D ---#
prim = zeros(5, 2)
prim[:, 1] .= [1.0, 1.0, 0.5, 3.0, 2.0]
prim[:, 2] .= [1 / 1864, 0.1, -0.3, 2.0, 1 / 1864]
# pure
KitCore.gauss_moments(prim[:, 1])
Mu1, Mv1, Mw1, MuL1, MuR1 = KitCore.gauss_moments(prim[:, 1], 0)
Mu2, Mv2, Mw2, MuL2, MuR2 = KitCore.gauss_moments(prim[:, 2], 0)
# mixture
Mu, Mv, Mw, MuL, MuR = KitCore.mixture_gauss_moments(prim, 0)

@test Mu[:, 1] ≈ Mu1 atol = 0.01
@test Mv[:, 1] ≈ Mv1 atol = 0.01
@test Mw[:, 1] ≈ Mw1 atol = 0.01
@test MuL[:, 1] ≈ MuL1 atol = 0.01
@test MuR[:, 1] ≈ MuR1 atol = 0.01

@test Mu[:, 2] ≈ Mu2 atol = 0.01
@test Mv[:, 2] ≈ Mv2 atol = 0.01
@test Mw[:, 2] ≈ Mw2 atol = 0.01
@test MuL[:, 2] ≈ MuL2 atol = 0.01
@test MuR[:, 2] ≈ MuR2 atol = 0.01

#--- generalized ---#
Mu1, Mv1, Mw1, MuL1, MuR1 = KitCore.gauss_moments(prim[:, 1], 0)
KitCore.moments_conserve_slope(zeros(5), Mu1, Mv1, Mw1, 2, 0, 0)

f1 = rand(16)
u1 = collect(-1:1/7.5:1)
w1 = ones(16)
KitCore.discrete_moments(f1, u1, w1, 1)
KitCore.moments_conserve(f1, u1, w1)
KitCore.moments_conserve(f1, f1, u1, w1)
KitCore.moments_conserve(f1, f1, f1, f1, u1, w1)

f2 = rand(16, 16)
u2 = randn(16, 16)
w2 = rand(16, 16)
KitCore.discrete_moments(f2, u2, w2, 1)
KitCore.moments_conserve(f2, u2, u2, w2)
KitCore.moments_conserve(f2, f2, u2, u2, w2)
KitCore.moments_conserve(f2, f2, f2, u2, u2, w2)

f3 = rand(16, 16, 16)
u3 = randn(16, 16, 16)
w3 = rand(16, 16, 16)
KitCore.discrete_moments(f3, u3, w3, 1)
KitCore.moments_conserve(f3, u3, u3, u3, w3)

KitCore.polyatomic_moments_conserve(f1, f1, f1, u1, w1)
KitCore.polyatomic_moments_conserve(f2, f2, f2, u2, u2, w2, KC.VDF{3,2})
KitCore.polyatomic_moments_conserve(f3, f3, u3, u3, u3, w3, KC.VDF{2,3})

KitCore.pressure(rand(3))
KitCore.pressure(f1, rand(3), u1, w1)
KitCore.pressure(f1, f1, rand(3), u1, w1, 2)
KitCore.pressure(f2, f2, rand(4), u2, u2, w2, 1, KC.VDF{2,2})
KitCore.pressure(f3, rand(5), u3, u3, u3, w3, 1, KC.VDF{1,3})
KitCore.stress(f1, rand(3), u1, w1)
KitCore.stress(f2, rand(4), u2, u2, w2)
KitCore.stress(f3, rand(5), u3, u3, u3, w3)

KitCore.heat_flux(f1, rand(3), u1, w1)
KitCore.heat_flux(f1, f1, rand(3), u1, w1)
KitCore.heat_flux(f1, f1, f1, rand(3), u1, w1)
KitCore.heat_flux(f2, rand(4), u2, u2, w2)
KitCore.heat_flux(f2, f2, rand(4), u2, u2, w2)
KitCore.heat_flux(f2, f2, f2, rand(4), u2, u2, w2)
KitCore.heat_flux(f3, rand(5), u3, u3, u3, w3)

# multi-species polyatomic
KitCore.mixture_polyatomic_moments_conserve(
    rand(8, 2),
    rand(8, 2),
    rand(8, 2),
    randn(8, 2),
    rand(8, 2),
)
KitCore.mixture_polyatomic_moments_conserve(
    rand(8, 2),
    rand(8, 2),
    rand(8, 2),
    randn(8, 2),
    randn(8, 2),
    rand(8, 2),
    KC.VDF{3,2},
)
KitCore.mixture_polyatomic_moments_conserve(
    rand(8, 2),
    rand(8, 2),
    randn(8, 2),
    randn(8, 2),
    randn(8, 2),
    rand(8, 2),
    KC.VDF{2,3},
)
