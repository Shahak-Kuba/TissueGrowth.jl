using TissueGrowth

# User input variables
D = 1.0
kf = 0.0005
A = 0.0
ρ₀ = 24.1796
Tmax = 26.0
r₀ = 1.05
btype = "hex"

# running simulation
θ,r,ρ,κ,σ,η = TissueGrowth.FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype);

# plotting solution
cmap = :jet
f = TissueGrowth.plotContinuumResults_Polar(θ, r, ρ, cmap, D, kf)