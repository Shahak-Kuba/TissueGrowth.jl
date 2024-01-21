# Discrete simulation functions and code description

# Index
```@index
```

# Data Structure
```@docs
TissueGrowth.SimResults_t
```

# 1D simulation code:
```@docs
TissueGrowth.sim1D()
```

# 2D simulation code:

```@docs
TissueGrowth.sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,prolif,death,embed,α,β,γ,event_δt,seed)
```

# ODE Problem:

```@docs
TissueGrowth.ODE_fnc_1D!(du,u,p,t) 
TissueGrowth.ODE_fnc_2D!(du,u,p,t) 
```

# Equations for Mechanical Relaxation and Normal Velocity:

```@docs
TissueGrowth.Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)
TissueGrowth.Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)
TissueGrowth.Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)
```

# Cell behaviour functions:

```@docs
TissueGrowth.P(event,ρ,α)
TissueGrowth.A(event,ρ,β)
TissueGrowth.E(event,ρ,γ)
TissueGrowth.affect!(integrator)
TissueGrowth.store_embed_cell_pos(pos)
```

# Custom modifier & Special functions
```@docs
Base.insert!
Base.deleteat!
TissueGrowth.lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ)
TissueGrowth.κ(rᵢ₋₁, rᵢ, rᵢ₊₁)
TissueGrowth.Ω(p)
```