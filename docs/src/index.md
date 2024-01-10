```@meta
CurrentModule = TissueGrowth
```

# TissueGrowth.jl

Documentation for [TissueGrowth]

```@index
```

# Data Structure
```@docs
TissueGrowth.SimResults_t
```

# 1D simulation code:
```@docs
TissueGrowth.SetupODEproblem1D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, dist_type)
TissueGrowth.sim1D()
```

# 2D simulation code:

```@docs
TissueGrowth.SetupODEproblem2D(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,growth_dir,prolif,death,embed,α,β,γ,dist_type)
TissueGrowth.sim2D()
```

# ODE Problem:

```@docs
TissueGrowth.ODE_fnc_1D!(du,u,p,t) 
TissueGrowth.ODE_fnc_2D!(du,u,p,t) 
```

# Equations for ODE problem:

```@docs
TissueGrowth.δ(rᵢ₊₁, rᵢ)
TissueGrowth.ρ(rᵢ₊₁, rᵢ)
TissueGrowth.τ(rᵢ₊₁, rᵢ₋₁)
TissueGrowth.n(rᵢ₊₁, rᵢ₋₁, type)
TissueGrowth.Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)
TissueGrowth.Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)
TissueGrowth.Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)
TissueGrowth.κ(rᵢ₋₁, rᵢ, rᵢ₊₁)
```

# Cell behaviour functions:

```@docs
TissueGrowth.P(event,ρ,α)
TissueGrowth.A(event,ρ,β)
TissueGrowth.E(event,ρ,γ)
TissueGrowth.affect!(integrator)
```


