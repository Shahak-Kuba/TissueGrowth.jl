module TissueGrowth
    # PACKAGES USED for solving equations
    using Base
    using DifferentialEquations
    using LinearAlgebra
    using Random
    using ElasticArrays
    using QuadGK
    using Roots
    # PACKAGES USED for benchmarking
    using BenchmarkTools
    # PACKAGES USED for plotting
    using Plots
    using Makie
    using CairoMakie
    using ColorSchemes
    using Colors
    # PACKAGES USED for misc
    using Printf
    using JLD2
    import FilePaths

    # DEVELOPED SIMULATION CODE

    # discrete simulation code
    include("Discrete/MechanicalEqns.jl")
    include("Discrete/CellBehaviours.jl")
    include("Discrete/PoreBoundariesV2.jl")
    include("Discrete/ProblemSetup.jl")
    include("Discrete/PlottingFncs1D.jl")
    include("Discrete/PlottingFncs2D.jl")
    include("Discrete/Simulation_1D.jl")
    include("Discrete/Simulation_2D.jl")
    include("Discrete/TissueGrowthODEproblem.jl")
    include("Discrete/GeometrySolvers.jl")
    include("Discrete/Misc.jl")
    include("Discrete/DataStructs.jl")
    include("Discrete/PostSimulation.jl")
    include("Discrete/ModifierFncs.jl")

    # continuum limit simulation code
    include("Continuum/Semi-Implicit_FD/FD_ContinuumSolvers.jl")
    include("Continuum/Semi-Implicit_FD/FD_SolverFncs.jl")
    include("Continuum/FVM_K-T/FVM_ContinuumSolver.jl")
    include("Continuum/FVM_K-T/FVM_SolverFncs.jl")
    include("Continuum/PlottingFncsPDE.jl")

    # including for comparion plotting 
    include("../Run/Comparison_Sims/ComparisonPlottingFncs.jl")
end
