using TissueGrowth
include("FD_ContinuumSolvers.jl")
include("../PlottingFncsPDE.jl")

# simulation parameters
D = 0.15;
kf = 0.01;
A = 0.00;
ρ₀ = 17.18247918322778;
growth_dir = "inward"
Tmax = 30.0
Xmax = 2π

x,h,ρ = FD_SolveContinuumLim_Cartesian(D,kf,A,ρ₀,growth_dir,Tmax,Xmax);

cmap = :jet
f1 = plotContinuumResults_Cartesian(x, h, ρ, D, kf, cmap)


# simulation parameters
D = 0.05;
kf = 0.001;
A = 0.00;
ρ₀ = 20;
growth_dir = "inward"
btype = "hex"
Tmax = 21
r₀ = 1.05

θ,R,ρ = FD_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype,growth_dir);

plotContinuumResults_Polar(θ, R, ρ, cmap, D, kf)

