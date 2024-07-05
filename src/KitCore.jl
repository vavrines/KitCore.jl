"""
KitCore.jl: The lightweight module of physical formulations in Kinetic.jl

Copyright (c) 2020-2024 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitCore

const KC = KitCore

import NonlinearSolve

using CUDA
using Distributions
using DocStringExtensions
using FFTW
using LinearAlgebra
using MultivariatePolynomials
using OffsetArrays
using Optim
using SpecialFunctions

using FastGaussQuadrature: gausslegendre
using Parameters: @with_kw
using Roots: Order1, find_zero
using SciMLNLSolve: NLSolveJL
using TypedPolynomials: Variable, @polyvar

include("Data/data.jl")
include("Macro/macro.jl")
include("Struct/struct.jl")
include("IO/io.jl")
include("Math/math.jl")
include("Geometry/geometry.jl")
include("Phase/phase.jl")
include("Theory/theory.jl")

export KC

end
