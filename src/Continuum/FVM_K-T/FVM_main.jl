include("../PlottingFncsPDE.jl")
include("FVM_ContinuumSolver.jl")

# User input variables
D = 0.015
kf = 0.001
A = 0.0
ρ₀ = 20.0
Tmax = 31
r₀ = 2.0

# running simulation
r,σ,η,ρ,θ = FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀);

# plotting solution
cmap = :jet
f = plotContinuumResults_Polar(θ, r, ρ, cmap)