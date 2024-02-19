# Example 2D continuum simulation 

## Cartesian simulation code:

%```@example sim2D_contiuum_cartesian
%using TissueGrowth
%# simulation parameters
%D = 0.0015;
%kf = 0.005;
%A = 0.00;
%ρ₀ = 14;
%growth_dir = "inward"
%Tmax = 50.0
%Xmax = 2π

%x,h,ρ = TissueGrowth.FD_SolveContinuumLim_Cartesian(D,kf,A,ρ₀,growth_dir,Tmax,Xmax);

%cmap = :jet
%f = TissueGrowth.plotContinuumResults_Cartesian(x, h, ρ, cmap, D, kf)
%```


## Polar simulation code:

```@example sim2D_contiuum_cartesian
using TissueGrowth
# simulation parameters
D = 0.0015;
kf = 0.001;
A = 0.00;
ρ₀ = 14;
growth_dir = "inward";
btype = "hex"
Tmax = 30.0;
r₀ = 1.05

θ,R,ρ = TissueGrowth.FD_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype,growth_dir);

cmap = :jet
f = TissueGrowth.plotContinuumResults_Polar(θ, R, ρ, cmap, D, kf)
```