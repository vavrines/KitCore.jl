mesh_quadrature(-5, -5, 8; type = "rectangle")
try
    mesh_quadrature(-5, -5, 8; type = "newton")
catch
    mesh_quadrature(-5, -5, 9; type = "newton")
end
mesh_quadrature(-5, -5, 8; type = "algebra")
mesh_quadrature(-5, -5, 8, -5, -5, 8; type = "rectangle")
mesh_quadrature(-5, -5, 9, -5, -5, 9; type = "newton")
mesh_quadrature(-5, -5, 8, -5, -5, 8; type = "algebra")
mesh_quadrature(-5, -5, 8, -5, -5, 8; type = "maxwell")
mesh_quadrature(-5, -5, 8, -5, -5, 8, -5, -5, 8; type = "rectangle")
mesh_quadrature(-5, -5, 9, -5, -5, 9, -5, -5, 9; type = "newton")
mesh_quadrature(-5, -5, 8, -5, -5, 8, -5, -5, 8; type = "algebra")

# spherical quadrature
KC.legendre_quadrature(6)
KC.octa_quadrature(8)
