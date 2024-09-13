module DiffAnalysis_AIC2024

using Reexport

# Asap dependencies
@reexport using Asap

# Analysis dependencies
using SparseArrays
@reexport using LinearSolve
@reexport using LinearAlgebra
@reexport using Statistics

# Optimization
using ChainRulesCore
@reexport using Zygote
@reexport using FiniteDifferences

# IO
@reexport using CSV, JSON, JLD2

# Data structures
include("Types/Types.jl")

# Utility functions
include("Utilities/Utilities.jl")

# Main functions
include("Functions/Functions.jl")

# Generators
include("Generators/Generators.jl")

# Visualization
@reexport using GLMakie
include("Visualization/Visualization.jl")

# Optimization
import Nonconvex, NonconvexNLopt
@reexport using Metaheuristics
include("Optimization/Optimization.jl")
export gradient_descent
export find_bands
export check_cstr
export constrained_optimization, unconstrained_optimization
export constrained_optimization_FD, unconstrained_optimization_FD

include("IO/GHAsap.jl")
export GHsave
end # module DiffAnalysis_AIC2024
