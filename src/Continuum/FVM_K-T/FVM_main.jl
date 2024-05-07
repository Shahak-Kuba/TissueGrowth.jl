using TissueGrowth
include("../PlottingFncsPDE.jl")
include("FVM_ContinuumSolver.jl")
include("FVM_SolverFncs.jl")

# User input variables
D = 0.001
kf = 20
A = 0.00
ρ₀ = 0.05
Tmax = 15
r₀ = 50
btype = "square"
growth_dir = "inward"

# running simulation
θ,r,ρ,κ,σ,η = FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype,growth_dir);

# plotting solution
cmap = :jet
f = plotContinuumResults_Polar(θ, r, ρ, cmap, D, kf)  