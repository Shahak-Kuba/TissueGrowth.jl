using TissueGrowth

# simulation parameters
D = 0.15;
kf = 0.01;
A = 0.00;
ρ₀ = 17.18247918322778;
growth_dir = "inward"
Tmax = 30.0
Xmax = 2π

x,h,ρ = TissueGrowth.SolveContinuumLim_Cartesian(D,kf,A,ρ₀,growth_dir,Tmax,Xmax);

cmap = :jet
f1 = TissueGrowth.plotContinuumResults_Cartesian(x, h, ρ, D, kf, cmap)


# simulation parameters
D = 0.015 ;
kf = 0.00053;
A = 0.00;
ρ₀ = 20;
growth_dir = "inward"
Tmax = 50.0

θ,R,ρ = TissueGrowth.SolveContinuumLim_Polar(D,kf,A,ρ₀,growth_dir,Tmax);

f2 = TissueGrowth.plotContinuumResults_Polar(θ, R, ρ, cmap)
