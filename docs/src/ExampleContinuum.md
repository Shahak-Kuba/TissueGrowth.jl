# Example 2D continuum simulation 

## Cartesian simulation code:

```@example sim2D_contiuum_cartesian
using TissueGrowth
# simulation parameters
D = 0.0015;
kf = 0.005;
A = 0.00;
ρ₀ = 14;
growth_dir = "inward"
Tmax = 50.0
Xmax = 2π

x,h,ρ = TissueGrowth.SolveContinuumLim_Cartesian(D,kf,A,ρ₀,growth_dir,Tmax,Xmax);

cmap = :jet
f = TissueGrowth.plotContinuumResults_Cartesian(x, h, ρ, cmap)
```


## Polar simulation code:

```@example sim2D_contiuum_cartesian
using TissueGrowth
# simulation parameters
D = 0.0015;
kf = 0.001;
A = 0.00;
ρ₀ = 14;
growth_dir = "inward";
Tmax = 30.0;

θ,R,ρ = TissueGrowth.SolveContinuumLim_Polar(D,kf,A,ρ₀,growth_dir,Tmax)

cmap = :jet
f = TissueGrowth.plotContinuumResults_Polar(θ, R, ρ, cmap)
```