# 2D
w = [1.0, 0.1, 0.3, 1.0]
cosa = cos(0.5)
sina = sin(0.5)
KC.global_frame(randn(2), cosa, sina)
KC.global_frame(randn(2), [cosa, sina])
KC.global_frame([1, 0], cosa, sina)
KC.local_frame(randn(2), cosa, sina)
KC.local_frame(randn(2), [cosa, sina])
KC.local_frame([1, 0], cosa, sina)
KC.global_frame(w, cosa, sina)
KC.local_frame(w, cosa, sina)

KC.global_frame(rand(3), rand(3, 3))
KC.local_frame(rand(3), rand(3, 3))
KC.global_frame(rand(5), rand(3, 3))
KC.local_frame(rand(5), rand(3, 3))

KC.local_velocity(rand(), rand(), [cosa, sina])
KC.local_velocity(rand(), rand(), [cosa, sina])

KC.unit_normal(rand(2), rand(2))
KC.unit_normal(rand(3), rand(3), rand(3))
KC.point_distance(rand(2), rand(2), rand(2))

x0 = 0.0
nx = 15
dx = 0.1
KC.uniform_mesh(x0, nx, dx)

x = collect(0:0.1:1)
y = collect(0:0.1:1)
z = collect(0:0.1:1)
KC.ndgrid(x)
KC.ndgrid(x, y)
KC.ndgrid(x, y, z)
KC.meshgrid(x)
KC.meshgrid(x, y)
KC.meshgrid(x, y, z)

KC.find_idx(randn(20), 0.13, mode = :uniform)
KC.find_idx(randn(20), 0.13, mode = :nonuniform)
