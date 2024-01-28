using TissueGrowth

# User input variables
D = 0.1
kf = 0.001
A = 0.0
ρ₀ = 20.0
Tmax = 31
r₀ = 2.0
btype = "square"

# running simulation
r,σ,η,ρ,θ = TissueGrowth.FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype);

# plotting solution
cmap = :jet
f = TissueGrowth.plotContinuumResults_Polar(θ, r, ρ, cmap, D, kf)