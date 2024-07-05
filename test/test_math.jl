x0 = 0
x1 = 1
nx = 5

KitCore.linspace(x0, x1, nx)
KitCore.heaviside(x0)
KitCore.fortsign(x0, x1)

KC.convergence_order(1e-2, 1e-3)
KC.L1_error(rand(3), rand(3), 1e-2)
KC.L2_error(rand(3), rand(3), 1e-2)
KC.L∞_error(rand(3), rand(3), 1e-2)

x = randn(16)
y = randn(16)
KitCore.@nametuple x y # NamedTuple constructor

KitCore.lgwt(12, -1, 1)
