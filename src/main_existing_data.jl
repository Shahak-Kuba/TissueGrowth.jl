using OrdinaryDiffEq
using CairoMakie
using ColorSchemes
using Colors
using ElasticArrays
using LinearAlgebra
using BenchmarkTools
using Printf
using CurveFit
using JLD2
import FilePaths

include("MechanicalEqns.jl")
include("PoreBoundaries.jl")
include("PlottingFncs1D.jl")
include("PlottingFncs2D.jl")
include("Simulation_1D.jl")
include("Simulation_2D.jl")
include("TissueGrowthODEproblem.jl")
include("GeometrySolvers.jl")
include("Misc.jl")
include("DataStructs.jl")
include("PostSimulation.jl")

# 1D simulations
sols1D = LoadData("test", "1D_simulations")
f = plotResults1D_Velocity(sols1D[1].u, sols1D[1].Vₙ)
f = plotResults1D_Velocity(sols1D[2].u, sols1D[2].Vₙ)
f = plotResults1D_Velocity(sols1D[3].u, sols1D[3].Vₙ)
f = plotResults1D_Velocity(sols1D[4].u, sols1D[4].Vₙ)


# 2D simulations
## Area Loss simulations
sols2d_δt = LoadData("AreaLoss_Simulations", "2D_Simulations")
plotAreaVStime_δt(sols2d_δt[1][1], sols2d_δt[1][2], sols2d_δt[1][3], sols2d_δt[1][4], "triangle")
plotAreaVStime_δt(sols2d_δt[2][1], sols2d_δt[2][2], sols2d_δt[2][3], sols2d_δt[2][4], "square")
plotAreaVStime_δt(sols2d_δt[3][1], sols2d_δt[3][2], sols2d_δt[3][3], sols2d_δt[3][4], "hex")
plotAreaVStime_δt(sols2d_δt[4][1], sols2d_δt[4][2], sols2d_δt[4][3], sols2d_δt[4][4], "star")
plotAreaVStime_δt(sols2d_δt[5][1], sols2d_δt[5][2], sols2d_δt[5][3], sols2d_δt[5][4], "cross")

## Diffusivity Simulations
sols2D_δt = LoadData("Diffusivity_Simulations", "2D_Simulations")
cmap = :viridis
f = plotResults2D(sols2D[1][1].u, sols2D[1][1].Density, cmap)
f = plotResults2D(sols2D[1][2].u, sols2D[1][2].Density, cmap)
f = plotResults2D(sols2D[1][3].u, sols2D[1][3].Density, cmap)
f = plotResults2D(sols2D[1][4].u, sols2D[1][4].Density, cmap)
f = plotResults2D(sols2D[1][5].u, sols2D[1][5].Density, cmap)
f = plotResults2D(sols2D[1][6].u, sols2D[1][6].Density, cmap)

f = plotResults2D(sols2D[2][1].u, sols2D[2][1].Density, cmap)
f = plotResults2D(sols2D[2][2].u, sols2D[2][2].Density, cmap)
f = plotResults2D(sols2D[2][3].u, sols2D[2][3].Density, cmap)
f = plotResults2D(sols2D[2][4].u, sols2D[2][4].Density, cmap)
f = plotResults2D(sols2D[2][5].u, sols2D[2][5].Density, cmap)
f = plotResults2D(sols2D[2][6].u, sols2D[2][6].Density, cmap)

f = plotResults2D(sols2D[3][1].u, sols2D[3][1].Density, cmap)
f = plotResults2D(sols2D[3][2].u, sols2D[3][2].Density, cmap)
f = plotResults2D(sols2D[3][3].u, sols2D[3][3].Density, cmap)
f = plotResults2D(sols2D[3][4].u, sols2D[3][4].Density, cmap)
f = plotResults2D(sols2D[3][5].u, sols2D[3][5].Density, cmap)
f = plotResults2D(sols2D[3][6].u, sols2D[3][6].Density, cmap)