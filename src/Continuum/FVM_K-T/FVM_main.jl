include("../PlottingFncsPDE.jl")
include("FVM_ContinuumSolver.jl")

# User input variables
D = 0.015
kf = 0.001
A = 0.00
ρ₀ = 20.0
Tmax = 50
r₀ = 2.0
btype = "circle"

# running simulation
r,σ,η,ρ,θ = FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype);

# plotting solution
cmap = :jet
f = plotContinuumResults_Polar(θ, r, ρ, cmap)  