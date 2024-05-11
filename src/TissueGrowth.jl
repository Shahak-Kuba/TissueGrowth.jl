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
    include("Discrete/GeneralEquations.jl")
    include("Discrete/DataStructs.jl")
    include("Discrete/ModifierFncs.jl")
    include("Discrete/Misc.jl")
    include("Discrete/PoreBoundaries.jl")

    include("Discrete/Model/CellMechanics.jl")
    include("Discrete/Model/CellBehaviours.jl")
    include("Discrete/Model/TissueSecretion.jl")
    include("Discrete/Model/AnalyticSolution.jl")
    
    include("Discrete/ProblemSetup.jl")
    include("Discrete/TissueGrowthODEproblem.jl")
    include("Discrete/PostSimulation.jl")
    
    include("Discrete/GrowthSimulation.jl")
    
    include("Discrete/Plotting/GeneralPlotting.jl")
    include("Discrete/Plotting/EmbeddedPlotting.jl")
    include("Discrete/Plotting/HistogramPlotting.jl")
    include("Discrete/Plotting/InterfaceAnimation.jl")
    include("Discrete/Plotting/PlottingFncs1D.jl")
    include("Discrete/Plotting/PlottingFncs2D.jl")

    # continuum limit simulation code
    include("Continuum/Semi-Implicit_FD/FD_ContinuumSolvers.jl")
    include("Continuum/Semi-Implicit_FD/FD_SolverFncs.jl")

    include("Continuum/FVM_K-T/FVM_ContinuumSolver.jl")
    include("Continuum/FVM_K-T/FVM_SolverFncs.jl")
    
    include("Continuum/PlottingFncsPDE.jl")

    # including for comparion plotting 
    include("../Simulations/Paper1_Simulations/C-Sim_Code/Pore_Simulations/Discrete_Continuum_Comparison_Sims/ComparisonPlottingFncs.jl")
end
