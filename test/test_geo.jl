w = [1.0, 0.1, 0.3, 1.0]
cosa = cos(0.5)
sina = sin(0.5)
KitCore.global_frame(randn(2), cosa, sina)
KitCore.local_frame(randn(2), cosa, sina)
KitCore.global_frame(w, cosa, sina)
KitCore.local_frame(w, cosa, sina)

KitCore.global_frame(rand(3), rand(3, 3))
KitCore.local_frame(rand(3), rand(3, 3))
KitCore.global_frame(rand(5), rand(3, 3))
KitCore.local_frame(rand(5), rand(3, 3))

KitCore.unit_normal(rand(2), rand(2))
KitCore.unit_normal(rand(3), rand(3), rand(3))
KitCore.point_distance(rand(2), rand(2), rand(2))

x0 = 0.0
nx = 15
dx = 0.1
KitCore.uniform_mesh(x0, nx, dx)

x = collect(0:0.1:1)
y = collect(0:0.1:1)
z = collect(0:0.1:1)
KitCore.ndgrid(x)
KitCore.ndgrid(x, y)
KitCore.ndgrid(x, y, z)
KitCore.meshgrid(x)
KitCore.meshgrid(x, y)
KitCore.meshgrid(x, y, z)

KitCore.find_idx(randn(20), 0.13, mode = :uniform)
KitCore.find_idx(randn(20), 0.13, mode = :nonuniform)
