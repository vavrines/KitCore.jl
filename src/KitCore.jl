"""
KitCore.jl: The lightweight module of physical formulations in Kinetic.jl

Copyright (c) 2020-2024 Tianbai Xiao <tianbaixiao@gmail.com>
"""

module KitCore

const KC = KitCore

import NonlinearSolve
import Roots: Order1, find_zero
import SciMLNLSolve: NLSolveJL
using CUDA
using DocStringExtensions
using FFTW
using LinearAlgebra
using OffsetArrays
using Parameters
using SpecialFunctions

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
