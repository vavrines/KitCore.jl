"""
KitCore.jl: The lightweight module of physical formulations in Kinetic.jl

Copyright (c) 2020-2024 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitCore

const KC = KitCore

import NonlinearSolve
import Roots: Order1, find_zero
import SciMLNLSolve: NLSolveJL
using Base.Threads: @threads
using CUDA
#using Distributions
using DocStringExtensions
#using FastGaussQuadrature
using FFTW
using LinearAlgebra
#using MultivariatePolynomials
using OffsetArrays
#using Optim
using Parameters
using SpecialFunctions
#using StaticArrays

include("Data/data.jl")
include("Macro/macro.jl")
include("Struct/struct.jl")
include("IO/io.jl")
include("Math/math.jl")
include("Geometry/geometry.jl")
include("Phase/phase.jl")
include("Theory/theory.jl")
include("Reconstruction/reconstruction.jl")

export KC

end
