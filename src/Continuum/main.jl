using TissueGrowth
include("ContinuumSolvers.jl")
include("PlottingFncsPDE.jl")

# simulation parameters
D = 1;
kf = 0.00053;
A = 0.00;
ρ₀ = 30;
growth_dir = "inward"
Tmax = 250.0
Xmax = 2π

x,h,ρ = SolveContinuumLim_Cartesian(D,kf,A,ρ₀,growth_dir,Tmax,Xmax);

cmap = :jet
f1 = plotContinuumResults_Cartesian(x, h, ρ, D, kf, cmap)


# simulation parameters
D = 1 ;
kf = 0.00053;
A = 0.00;
ρ₀ = 30;
growth_dir = "inward"
Tmax = 250.0

θ,R,ρ = SolveContinuumLim_Polar(D,kf,A,ρ₀,growth_dir,Tmax);

plotContinuumResults_Polar(θ, R, ρ, cmap)
