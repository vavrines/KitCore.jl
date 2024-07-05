using Test, KitCore
using KitCore.SpecialFunctions

cd(@__DIR__)
include("test_data.jl")
include("test_macro.jl")
include("test_struct.jl")
include("test_io.jl")
include("test_math.jl")
include("test_geo.jl")
include("test_phase.jl")
include("test_theory.jl")
include("test_moments.jl")
include("test_reconstruction.jl")
