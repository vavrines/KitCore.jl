#--- continuum ---#
prim = [1.0, 0.0, 1.0]
KitCore.prim_conserve(prim, 3.0)
KitCore.prim_conserve(prim[1], prim[2], prim[3], 3.0)

mprim = hcat(prim, prim)
KitCore.mixture_prim_conserve(mprim, 3.0)

KitCore.conserve_prim(prim[1])
KitCore.conserve_prim(prim[1], 1.0)
KitCore.conserve_prim(prim, 3.0)
KitCore.conserve_prim(prim[1], prim[2], prim[3], 3.0)
KitCore.mixture_conserve_prim(mprim, 3.0)

# polyatomic
KitCore.prim_conserve([1.0, 0.0, 1.0, 1.0, 1.0], 5 / 3, 2)
KitCore.prim_conserve([1.0, 0.0, 0.0, 1.0, 1.0, 1.0], 5 / 3, 2)
KitCore.conserve_prim([1.0, 0.0, 1.0, 0.1], 5 / 3, 2)
KitCore.conserve_prim([1.0, 0.0, 0.0, 1.0, 0.1], 5 / 3, 2)
# multi-species polyatomic
KitCore.mixture_prim_conserve(rand(5, 2), 5 / 3, 2)
KitCore.mixture_prim_conserve(rand(6, 2), 5 / 3, 2)
KitCore.mixture_prim_conserve(rand(7, 2), 5 / 3, 2)
KitCore.mixture_conserve_prim(rand(4, 2), 2, 2)
KitCore.mixture_conserve_prim(rand(5, 2), 1, 2)
KitCore.mixture_conserve_prim(rand(6, 2), 0, 2)

prim = [1.0, 0.2, 0.3, -0.1, 1.0]
mprim = hcat(prim, prim)
KitCore.em_coefficients(mprim, randn(3), randn(3), 100, 0.01, 0.01, 0.001)

KitCore.advection_flux(1.0, -0.1)
KitCore.burgers_flux(1.0)
KitCore.euler_flux(prim, 3.0)
KitCore.euler_jacobi(prim, 3.0)

#--- thermo ---#
KitCore.heat_capacity_ratio(2.0, 1)
KitCore.heat_capacity_ratio(2.0, 2)
KitCore.heat_capacity_ratio(2.0, 3)
KitCore.heat_capacity_ratio(2.0, 2, 1)
KitCore.heat_capacity_ratio(2.0, 2, 2)
KitCore.heat_capacity_ratio(2.0, 2, 3)

KitCore.sound_speed(1.0, 5 / 3)
KitCore.sound_speed([1.0, 0.0, 1.0], 5 / 3)
KitCore.sound_speed(rand(3, 2), 5 / 3)

#--- atom ---#
KitCore.pdf_slope(1.0, 0.1)
KitCore.pdf_slope(prim, randn(5), 0.0)
KitCore.mixture_pdf_slope(mprim, randn(5, 2), 0.0)

u = collect(-5:0.2:5)
ω = ones(51) .* 0.2
prim = [1.0, 0.0, 1.0]

M = KitCore.maxwellian(u, prim)
KitCore.maxwellian(randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KitCore.maxwellian(
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    [1.0, 0.0, 0.0, 0.0, 1.0],
)

KitCore.energy_maxwellian(M, prim, 2)
KitCore.mixture_energy_maxwellian(hcat(M, M), hcat(prim, prim), 2)
KitCore.mixture_energy_maxwellian(rand(16, 16, 2), hcat(prim, prim), 2)

KitCore.maxwellian!(rand(16), rand(16), prim)
KitCore.maxwellian!(rand(16, 16), randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KitCore.maxwellian!(
    rand(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    [1.0, 0.0, 0.0, 0.0, 1.0],
)

mprim = hcat(prim, prim)
KitCore.mixture_maxwellian(hcat(u, u), mprim)
KitCore.mixture_maxwellian(randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KitCore.mixture_maxwellian(
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    rand(5, 2),
)

KitCore.mixture_maxwellian!(randn(8, 2), randn(8, 2), rand(3, 2))
KitCore.mixture_maxwellian!(randn(8, 8, 2), randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KitCore.mixture_maxwellian!(
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    rand(5, 2),
)

KitCore.shakhov(u, M, 0.01, prim, 1.0)
KitCore.shakhov(u, M, M, 0.01, prim, 1.0, 2.0)
KitCore.shakhov(randn(16, 16), randn(16, 16), rand(16, 16), rand(2), [1.0, 0.0, 1.0], 1.0)
KitCore.shakhov(
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 0.0, 1.0],
    1.0,
    1.0,
)
KitCore.shakhov(
    randn(8, 8, 8),
    randn(8, 8, 8),
    rand(8, 8, 8),
    rand(8, 8, 8),
    rand(3),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    1.0,
)

KitCore.shakhov!(randn(16), randn(16), rand(16), 0.01, prim, 1.0)
KitCore.shakhov!(randn(16), randn(16), randn(16), rand(16), rand(16), 0.01, prim, 1.0, 2.0)
KitCore.shakhov!(
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 1.0],
    1.0,
)
KitCore.shakhov!(
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 0.0, 1.0],
    1.0,
    1.0,
)
KitCore.shakhov!(
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    rand(8, 8, 8),
    rand(8, 8, 8),
    rand(3),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    1.0,
)

# ES-BGK
#vs = VSpace1D()
u, _, weights = mesh_quadrature(-5, 5, 16)
prim = [1.0, 0, 1]
f = maxwellian(u, prim)
KitCore.esbgk_ode!(zero(f), f, (u, weights, prim, 2 / 3, 1), 0.0)

#vs = VSpace2D()
u, v, _, _, weights = mesh_quadrature(-5, 5, 16, -5, 5, 16)
prim = [1.0, 0, 0, 1]
f = maxwellian(u, v, prim)
KitCore.esbgk_ode!(zero(f), f, (u, v, weights, prim, 2 / 3, 1), 0.0)

#vs = VSpace3D()
u, v, w, _, _, _, weights = mesh_quadrature(-5, 5, 16, -5, 5, 16, -5, 5, 16)
prim = [1.0, 0, 0, 0, 1]
f = maxwellian(u, v, w, prim)
KitCore.esbgk_ode!(zero(f), f, (u, v, w, weights, prim, 2 / 3, 1), 0.0)

# Rykov
KitCore.polyatomic_maxwellian!(
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    rand(5),
    4,
    2,
)
KitCore.polyatomic_maxwellian!(
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    randn(8, 8),
    randn(8, 8),
    rand(6),
    4,
    2,
)

KitCore.rykov!(
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(2),
    rand(5),
    0.72,
    4,
    1 / 1.55,
    0.2354,
    0.3049,
)
KitCore.rykov!(
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    randn(8, 8),
    randn(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(4),
    rand(6),
    0.72,
    4,
    1 / 1.55,
    0.2354,
    0.3049,
)

# BIP
KitCore.polyatomic_maxwellian!(zeros(16), zeros(16), zeros(16), randn(16), rand(5), 2, 2)
KitCore.polyatomic_maxwellian!(
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    randn(16),
    rand(6),
    1,
    2,
)
KitCore.polyatomic_maxwellian!(
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    randn(16),
    randn(16),
    rand(7),
    1,
    2,
)

# multi-species BIP
KitCore.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    rand(5, 2),
    2,
    2,
)
KitCore.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    randn(16, 2),
    rand(6, 2),
    1,
    2,
)
KitCore.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    randn(16, 2),
    randn(16, 2),
    rand(7, 2),
    0,
    2,
)

KitCore.f_maxwellian(rand(16), randn(16), rand(16))
KitCore.f_maxwellian(rand(16), rand(16), randn(16), rand(16))

KitCore.reduce_distribution(randn(16, 51), ω, 1)
KitCore.reduce_distribution(randn(16, 24, 24), ones(24, 24), 1)
KitCore.reduce_distribution(
    randn(16, 24, 24),
    randn(16, 24, 24),
    randn(16, 24, 24),
    ones(24, 24),
    1,
)
KitCore.full_distribution(M, M, randn(length(ω)), ω, ones(51, 24, 24), ones(51, 24, 24), 1.0, 3.0)
KitCore.full_distribution(M, M, randn(length(ω)), ω, ones(51, 24, 24), ones(51, 24, 24), prim, 3.0)

KitCore.ref_vhs_vis(1.0, 1.0, 0.5)
KitCore.vhs_collision_time(prim, 1e-3, 0.81)

KitCore.νbgk_relaxation_time(0.1, 1, rand(3), KC.Class{1})
KitCore.νbgk_relaxation_time(0.1, 1, rand(3), KC.Class{2})
KitCore.νbgk_relaxation_time(0.1, 1, 1, rand(4), KC.Class{3})
KitCore.νbgk_relaxation_time(0.1, 1, 1, 1, rand(5), KC.Class{4})

KitCore.νshakhov_relaxation_time(0.1, 1, rand(3))

KitCore.rykov_zr(100, 91.5, 18.1)

KitCore.hs_boltz_kn(1e-3, 1.0)
#vs = VSpace3D(-5, 5, 16, -5, 5, 16, -5, 5, 16)
#fsm = KitCore.fsm_kernel(vs, 1e-3)
fsm = kernel_mode(5, 8, 8, 8, 48, 28, 28, 1)
phi, psi, phipsi =
    KitCore.kernel_mode(5, 5.0, 5.0, 5.0, 0.1, 0.1, 0.1, 16, 16, 16, 1.0, quad_num = 16)
KitCore.kernel_mode(5, 5.0, 5.0, 5.0, 16, 16, 16, 1.0, quad_num = 16)
KitCore.kernel_mode(5, 5.0, 5.0, 0.1, 0.1, 16, 16, quad_num = 16)
#KitCore.boltzmann_fft(rand(16, 16, 16), fsm)
#KitCore.boltzmann_fft!(rand(16, 16, 16), rand(16, 16, 16), fsm)

KitCore.boltzmann_ode!(zeros(16, 16, 16), rand(16, 16, 16), (1.0, 5, phi, psi, phipsi), 0.0)
KitCore.bgk_ode!(zeros(16, 16, 16), rand(16, 16, 16), (rand(16, 16, 16), 1e-2), 0.0)

#=vs = KitCore.VSpace3D(-5.0, 5.0, 16, -5.0, 5.0, 16, -5.0, 5.0, 16, type = "algebra")
u, v, w = vs.u[:, 1, 1], vs.v[1, :, 1], vs.w[1, 1, :]
vnu = hcat(vs.u[:], vs.v[:], vs.w[:])
uuni1d = linspace(vs.u[1, 1, 1], vs.u[end, 1, 1], 16)
vuni1d = linspace(vs.v[1, 1, 1], vs.v[1, end, 1], 16)
wuni1d = linspace(vs.w[1, 1, 1], vs.w[1, 1, end], 16)
u13d = [uuni1d[i] for i = 1:16, j = 1:16, k = 1:16]
v13d = [vuni1d[j] for i = 1:16, j = 1:16, k = 1:16]
w13d = [wuni1d[k] for i = 1:16, j = 1:16, k = 1:16]
vuni = hcat(u13d[:], v13d[:], w13d[:])=#

τ = KitCore.aap_hs_collision_time(mprim, 1.0, 0.5, 0.5, 0.5, 1.0)
KitCore.aap_hs_prim(mprim, τ, 1.0, 0.5, 0.5, 0.5, 1.0)
KitCore.aap_hs_prim(rand(4, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)
KitCore.aap_hs_prim(rand(5, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)

KitCore.aap_hs_diffeq!(
    similar(mprim),
    mprim,
    [τ[1], τ[2], 1.0, 0.5, 0.5, 0.5, 1.0, 3.0],
    0.0,
)
KitCore.shift_pdf!(M, 1.0, 1e-4, 1e-4)
KitCore.shift_pdf!(rand(16, 2), randn(2), rand(2), 1e-4)

KitCore.chapman_enskog(rand(16), [1.0, 0.0, 1.0], rand(3), rand(3), 0.1)
KitCore.chapman_enskog(rand(16), [1.0, 0.0, 1.0], zeros(3), 0, 0.1)
KitCore.chapman_enskog(
    rand(16, 16),
    rand(16, 16),
    [1.0, 0.0, 0.0, 1.0],
    rand(4),
    rand(4),
    rand(4),
    0.1,
)
KitCore.chapman_enskog(
    rand(16, 16),
    rand(16, 16),
    [1.0, 0.0, 0.0, 1.0],
    zeros(4),
    zeros(4),
    0.0,
    0.1,
)
KitCore.chapman_enskog(
    rand(16, 16, 16),
    rand(16, 16, 16),
    rand(16, 16, 16),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    rand(5),
    rand(5),
    rand(5),
    rand(5),
    0.1,
)
KitCore.chapman_enskog(
    rand(16, 16, 16),
    rand(16, 16, 16),
    rand(16, 16, 16),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    zeros(5),
    zeros(5),
    zeros(5),
    0.0,
    0.1,
)

#=vs = KitCore.VSpace1D()
KitCore.collision_invariant(rand(3), vs)
KitCore.collision_invariant(rand(3, 2), vs)
vs2 = KitCore.VSpace2D()
KitCore.collision_invariant(rand(4), vs2)
vs3 = KitCore.VSpace3D()
KitCore.collision_invariant(rand(5), vs3)=#

#--- quantum ---#
u, _, weights = mesh_quadrature(-5, 5, 16)
f0 = 0.5 * (1 / π)^0.5 .* (exp.(-(u .- 0.99) .^ 2) .+ exp.(-(u .+ 0.99) .^ 2))
w0 = moments_conserve(f0, u, weights)
prim0 = quantum_conserve_prim(w0, 2, :fd)
quantum_prim_conserve(prim0, 2, :fd)
prim0 = quantum_conserve_prim(w0, 2, :be)
quantum_prim_conserve(prim0, 2, :be)

u, v, _, _, weights = mesh_quadrature(-5, 5, 16, -5, 5, 16)
f0 =
    0.5 * (1 / π)^0.5 .* (exp.(-(u .- 0.99) .^ 2) .+ exp.(-(u .+ 0.99) .^ 2)) .*
    exp.(-v .^ 2)
w0 = moments_conserve(f0, u, v, weights)
prim0 = quantum_conserve_prim(w0, 2, :fd)
quantum_prim_conserve(prim0, 2, :fd)
prim0 = quantum_conserve_prim(w0, 2, :be)
quantum_prim_conserve(prim0, 2, :be)

#--- Hermite ---#
#vs = VSpace1D(-5, 5, 36)
u, _, weights = mesh_quadrature(-5, 5, 16)
prim = [2.0, 0.5, 0.6]
f = maxwellian(u, prim)
df = hermite_force(f, u, weights, prim, 11, 1.0)

#--- Riemann solution ---#
KC.sample_riemann_solution(
    [-0.5, -0.2, 0.1, 0.3, 0.5],
    0.2,
    KC.HydroStatus(1.0, 0.0, 1.0, 1.4),
    KC.HydroStatus(0.125, 0.0, 0.1, 1.4),
)
