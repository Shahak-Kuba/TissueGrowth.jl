```@meta
CurrentModule = TissueGrowth
```

# TissueGrowth.jl

TissueGrowth.jl is a open source project package developed to simulate the evoluition of biological tissue interface during tissue growth. Our model considers mechanical interactions between neighbouring cells and a secretion rate of new tissue material that is proportional to cell density. In this project we include cell proliferation, apoptosis (death) and embedment as stochastic processes. Below is a discription of our model.

# Model Description

This cell-based mathematical model 

Calculates the mechanical relaxation and normal velocity in a system with periodic boundary conditions. The derivatives are based on the discrete equation given by:

```math
\frac{\text{d}\mathbf{u}ᵢ}{\text{d}t} = \frac{1}{η}\bigg((\mathbf{F}ₛ⁺ - \mathbf{F}ₛ⁻)\cdot\bm{τ}\bigg)\bm{τ} + Vₙ\mathbf{n} 
```
where,

```math
\mathbf{F}ₛ⁺ = f(\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert) \frac{\mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ}{\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert}, \hspace{0.5cm} \mathbf{F}ₛ⁻ = f(\lVert \mathbf{u}ᵢ - \mathbf{u}ᵢ₋₁ \rVert) \frac{\mathbf{u}ᵢ - \mathbf{u}ᵢ₋₁}{\lVert \mathbf{u}ᵢ - \mathbf{u}ᵢ₋₁ \rVert}
```
In this current version of the simulation code we use a nonlinear resorting force given by,

```math
f(\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert) = kₛ \bigg(\frac{1}{l_{0}} - \frac{1}{\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert }\bigg)
```
and a normal velocity which is proportional to density such that,
```math
Vₙ = k_{f}\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert
```

given $f(\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert)$ is the restoration force function and $Vₙ$ is the velocity in the normal direction.

# Documentation Index
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

# Custom modifier functions
```@docs
Base.insert!
Base.deleteat!
```
