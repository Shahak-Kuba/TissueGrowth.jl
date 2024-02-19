# Contiunuum model (PDE solver)

For mid to high diffusivities ($0.05\leq D \leq 1$) we use the Finite Difference Upwinding (in space) scheme, and for low to mid diffusivities ($0\leq D \leq 0.05$) we use FVM following the Kruganov-Tadmor scheme (Alias & Buenzli, 2017). Note: The FVM solver is only implemented for solving domains in polar coordinates.

# Finite Difference Solver Functions
## PDE in Cartesian coordinates
```@docs
TissueGrowth.FD_SolveContinuumLim_Cartesian(D, kf, A, ρ₀, growth_dir, Tmax, Xmax)
```

## PDE in Polar coordinates
```@docs
TissueGrowth.FD_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype,growth_dir)
```

# FVM Kruganov-Tadmor Solver Functions
```@docs
TissueGrowth.FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax,r₀,btype, growth_dir)
```