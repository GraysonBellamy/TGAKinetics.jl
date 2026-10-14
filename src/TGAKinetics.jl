"""
    TGAKinetics

Fit solid-phase reaction schemes to nitrogen TGA data at one or more heating
rates. Everything is set in one TOML file: the data files, how to read them,
the reaction scheme, which parameters are free (and their bounds), the
objective weights and the optimizer settings.

Entry points: [`check`](@ref), [`fit`](@ref), [`simulate`](@ref). The
command-line wrapper is `bin/tgafit.jl`.
"""
module TGAKinetics

using ADTypes: AutoFiniteDiff
using CSV
using Dates
using DiffEqBase: DEVerbosity
using LinearAlgebra
using OrdinaryDiffEqRosenbrock: Rodas5P, Rosenbrock23
using Printf
using Random
using SciMLBase: ODEFunction, ODEProblem, init, solve!, successful_retcode
using SciMLLogging: SciMLLogging
using Sobol
using SparseArrays
using StaticArrays
using Statistics: mean
using TOML

const R_GAS = 8.31446261815324  # J/(mol K)

include("config.jl")
include("scheme.jl")
include("smoothing.jl")
include("data.jl")
include("model.jl")
include("objective.jl")
include("de1220.jl")
include("output.jl")
include("run.jl")

export load_config, check, fit, simulate

end
