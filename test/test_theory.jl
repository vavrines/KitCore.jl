#--- closure.jl ---#
KC.moment_basis([1.0], 2)
KC.moment_basis([1.0], [1.0], 2)
KC.moment_basis([1.0], [1.0], [1.0], 2)

quadratureorder = 2
points, weights = KC.octa_quadrature(quadratureorder)
nq = size(points, 1)
L = 1
ne = (L + 1)^2
α = zeros(ne)
u = [2.0, 0.0, 0.0, 0.0]
m = KC.eval_spherharmonic(points, L)

KC.eval_sphermonomial(rand(6), L)
KC.eval_sphermonomial(points, L)

res = KC.optimize_closure(α, m, weights, u, KC.maxwell_boltzmann_dual)
u1 = KC.realizable_reconstruct(
    res.minimizer,
    m,
    weights,
    KC.maxwell_boltzmann_dual_prime,
)

#--- continuum ---#
prim = [1.0, 0.0, 1.0]
KC.prim_conserve(prim, 3.0)
KC.prim_conserve(prim[1], prim[2], prim[3], 3.0)

mprim = hcat(prim, prim)
KC.mixture_prim_conserve(mprim, 3.0)

KC.conserve_prim(prim[1])
KC.conserve_prim(prim[1], 1.0)
KC.conserve_prim(prim, 3.0)
KC.conserve_prim(prim[1], prim[2], prim[3], 3.0)
KC.mixture_conserve_prim(mprim, 3.0)

# polyatomic
KC.prim_conserve([1.0, 0.0, 1.0, 1.0, 1.0], 5 / 3, 2)
KC.prim_conserve([1.0, 0.0, 0.0, 1.0, 1.0, 1.0], 5 / 3, 2)
KC.conserve_prim([1.0, 0.0, 1.0, 0.1], 5 / 3, 2)
KC.conserve_prim([1.0, 0.0, 0.0, 1.0, 0.1], 5 / 3, 2)
# multi-species polyatomic
KC.mixture_prim_conserve(rand(5, 2), 5 / 3, 2)
KC.mixture_prim_conserve(rand(6, 2), 5 / 3, 2)
KC.mixture_prim_conserve(rand(7, 2), 5 / 3, 2)
KC.mixture_conserve_prim(rand(4, 2), 2, 2)
KC.mixture_conserve_prim(rand(5, 2), 1, 2)
KC.mixture_conserve_prim(rand(6, 2), 0, 2)

prim = [1.0, 0.2, 0.3, -0.1, 1.0]
mprim = hcat(prim, prim)
KC.em_coefficients(mprim, randn(3), randn(3), 100, 0.01, 0.01, 0.001)

KC.advection_flux(1.0, -0.1)
KC.burgers_flux(1.0)
KC.euler_flux(rand(3), 3.0)
KC.euler_flux(rand(4), 3.0)
KC.euler_flux(prim, 3.0)
KC.euler_jacobi(prim, 3.0)

#--- thermo ---#
KC.heat_capacity_ratio(2.0, 1)
KC.heat_capacity_ratio(2.0, 2)
KC.heat_capacity_ratio(2.0, 3)
KC.heat_capacity_ratio(2.0, 2, 1)
KC.heat_capacity_ratio(2.0, 2, 2)
KC.heat_capacity_ratio(2.0, 2, 3)

KC.sound_speed(1.0, 5 / 3)
KC.sound_speed([1.0, 0.0, 1.0], 5 / 3)
KC.sound_speed(rand(3, 2), 5 / 3)

#--- atom ---#
KC.pdf_slope(1.0, 0.1)
KC.pdf_slope(prim, randn(5), 0.0)
KC.mixture_pdf_slope(mprim, randn(5, 2), 0.0)

u = collect(-5:0.2:5)
ω = ones(51) .* 0.2
prim = [1.0, 0.0, 1.0]

M = KC.maxwellian(u, prim)
KC.maxwellian(randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KC.maxwellian(
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    [1.0, 0.0, 0.0, 0.0, 1.0],
)

KC.energy_maxwellian(M, prim, 2)
KC.mixture_energy_maxwellian(hcat(M, M), hcat(prim, prim), 2)
KC.mixture_energy_maxwellian(rand(16, 16, 2), hcat(prim, prim), 2)

KC.maxwellian!(rand(16), rand(16), prim)
KC.maxwellian!(rand(16, 16), randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KC.maxwellian!(
    rand(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    [1.0, 0.0, 0.0, 0.0, 1.0],
)

mprim = hcat(prim, prim)
KC.mixture_maxwellian(hcat(u, u), mprim)
KC.mixture_maxwellian(randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KC.mixture_maxwellian(
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    rand(5, 2),
)

KC.mixture_maxwellian!(randn(8, 2), randn(8, 2), rand(3, 2))
KC.mixture_maxwellian!(randn(8, 8, 2), randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KC.mixture_maxwellian!(
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    rand(5, 2),
)

KC.shakhov(u, M, 0.01, prim, 1.0)
KC.shakhov(u, M, M, 0.01, prim, 1.0, 2.0)
KC.shakhov(randn(16, 16), randn(16, 16), rand(16, 16), rand(2), [1.0, 0.0, 1.0], 1.0)
KC.shakhov(
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 0.0, 1.0],
    1.0,
    1.0,
)
KC.shakhov(
    randn(8, 8, 8),
    randn(8, 8, 8),
    rand(8, 8, 8),
    rand(8, 8, 8),
    rand(3),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    1.0,
)

KC.shakhov!(randn(16), randn(16), rand(16), 0.01, prim, 1.0)
KC.shakhov!(randn(16), randn(16), randn(16), rand(16), rand(16), 0.01, prim, 1.0, 2.0)
KC.shakhov!(
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 1.0],
    1.0,
)
KC.shakhov!(
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
KC.shakhov!(
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
KC.esbgk_ode!(zero(f), f, (u, weights, prim, 2 / 3, 1), 0.0)

#vs = VSpace2D()
u, v, _, _, weights = mesh_quadrature(-5, 5, 16, -5, 5, 16)
prim = [1.0, 0, 0, 1]
f = maxwellian(u, v, prim)
KC.esbgk_ode!(zero(f), f, (u, v, weights, prim, 2 / 3, 1), 0.0)

#vs = VSpace3D()
u, v, w, _, _, _, weights = mesh_quadrature(-5, 5, 16, -5, 5, 16, -5, 5, 16)
prim = [1.0, 0, 0, 0, 1]
f = maxwellian(u, v, w, prim)
KC.esbgk_ode!(zero(f), f, (u, v, w, weights, prim, 2 / 3, 1), 0.0)

# Rykov
KC.polyatomic_maxwellian!(
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
KC.polyatomic_maxwellian!(
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

KC.rykov!(
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
KC.rykov!(
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
KC.polyatomic_maxwellian!(zeros(16), zeros(16), zeros(16), randn(16), rand(5), 2, 2)
KC.polyatomic_maxwellian!(
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    randn(16),
    rand(6),
    1,
    2,
)
KC.polyatomic_maxwellian!(
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
KC.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    rand(5, 2),
    2,
    2,
)
KC.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    randn(16, 2),
    rand(6, 2),
    1,
    2,
)
KC.mixture_polyatomic_maxwellian!(
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

KC.f_maxwellian(rand(16), randn(16), rand(16))
KC.f_maxwellian(rand(16), rand(16), randn(16), rand(16))

KC.reduce_distribution(randn(16, 51), ω, 1)
KC.reduce_distribution(randn(16, 24, 24), ones(24, 24), 1)
KC.reduce_distribution(
    randn(16, 24, 24),
    randn(16, 24, 24),
    randn(16, 24, 24),
    ones(24, 24),
    1,
)
KC.full_distribution(
    M,
    M,
    randn(length(ω)),
    ω,
    ones(51, 24, 24),
    ones(51, 24, 24),
    1.0,
    3.0,
)
KC.full_distribution(
    M,
    M,
    randn(length(ω)),
    ω,
    ones(51, 24, 24),
    ones(51, 24, 24),
    prim,
    3.0,
)

KC.ref_vhs_vis(1.0, 1.0, 0.5)
KC.vhs_collision_time(prim, 1e-3, 0.81)

KC.νbgk_relaxation_time(0.1, 1, rand(3), KC.Class{1})
KC.νbgk_relaxation_time(0.1, 1, rand(3), KC.Class{2})
KC.νbgk_relaxation_time(0.1, 1, 1, rand(4), KC.Class{3})
KC.νbgk_relaxation_time(0.1, 1, 1, 1, rand(5), KC.Class{4})

KC.νshakhov_relaxation_time(0.1, 1, rand(3))

KC.rykov_zr(100, 91.5, 18.1)

KC.hs_boltz_kn(1e-3, 1.0)
#vs = VSpace3D(-5, 5, 16, -5, 5, 16, -5, 5, 16)
#fsm = KC.fsm_kernel(vs, 1e-3)
fsm = kernel_mode(5, 8, 8, 8, 48, 28, 28, 1)
phi, psi, phipsi =
    KC.kernel_mode(5, 5.0, 5.0, 5.0, 0.1, 0.1, 0.1, 16, 16, 16, 1.0, quad_num = 16)
KC.kernel_mode(5, 5.0, 5.0, 5.0, 16, 16, 16, 1.0, quad_num = 16)
KC.kernel_mode(5, 5.0, 5.0, 0.1, 0.1, 16, 16, quad_num = 16)
#KC.boltzmann_fft(rand(16, 16, 16), fsm)
#KC.boltzmann_fft!(rand(16, 16, 16), rand(16, 16, 16), fsm)

KC.boltzmann_ode!(zeros(16, 16, 16), rand(16, 16, 16), (1.0, 5, phi, psi, phipsi), 0.0)
KC.bgk_ode!(zeros(16, 16, 16), rand(16, 16, 16), (rand(16, 16, 16), 1e-2), 0.0)

#=vs = KC.VSpace3D(-5.0, 5.0, 16, -5.0, 5.0, 16, -5.0, 5.0, 16, type = "algebra")
u, v, w = vs.u[:, 1, 1], vs.v[1, :, 1], vs.w[1, 1, :]
vnu = hcat(vs.u[:], vs.v[:], vs.w[:])
uuni1d = linspace(vs.u[1, 1, 1], vs.u[end, 1, 1], 16)
vuni1d = linspace(vs.v[1, 1, 1], vs.v[1, end, 1], 16)
wuni1d = linspace(vs.w[1, 1, 1], vs.w[1, 1, end], 16)
u13d = [uuni1d[i] for i = 1:16, j = 1:16, k = 1:16]
v13d = [vuni1d[j] for i = 1:16, j = 1:16, k = 1:16]
w13d = [wuni1d[k] for i = 1:16, j = 1:16, k = 1:16]
vuni = hcat(u13d[:], v13d[:], w13d[:])=#

τ = KC.aap_hs_collision_time(mprim, 1.0, 0.5, 0.5, 0.5, 1.0)
KC.aap_hs_prim(mprim, τ, 1.0, 0.5, 0.5, 0.5, 1.0)
KC.aap_hs_prim(rand(4, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)
KC.aap_hs_prim(rand(5, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)

KC.aap_hs_diffeq!(
    similar(mprim),
    mprim,
    [τ[1], τ[2], 1.0, 0.5, 0.5, 0.5, 1.0, 3.0],
    0.0,
)
KC.shift_pdf!(M, 1.0, 1e-4, 1e-4)
KC.shift_pdf!(rand(16, 2), randn(2), rand(2), 1e-4)

KC.chapman_enskog(rand(16), [1.0, 0.0, 1.0], rand(3), rand(3), 0.1)
KC.chapman_enskog(rand(16), [1.0, 0.0, 1.0], zeros(3), 0, 0.1)
KC.chapman_enskog(
    rand(16, 16),
    rand(16, 16),
    [1.0, 0.0, 0.0, 1.0],
    rand(4),
    rand(4),
    rand(4),
    0.1,
)
KC.chapman_enskog(
    rand(16, 16),
    rand(16, 16),
    [1.0, 0.0, 0.0, 1.0],
    zeros(4),
    zeros(4),
    0.0,
    0.1,
)
KC.chapman_enskog(
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
KC.chapman_enskog(
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

#=vs = KC.VSpace1D()
KC.collision_invariant(rand(3), vs)
KC.collision_invariant(rand(3, 2), vs)
vs2 = KC.VSpace2D()
KC.collision_invariant(rand(4), vs2)
vs3 = KC.VSpace3D()
KC.collision_invariant(rand(5), vs3)=#

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
