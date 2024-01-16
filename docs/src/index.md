```@meta
CurrentModule = TissueGrowth
```

# TissueGrowth.jl

TissueGrowth.jl is a open source project package developed to simulate the evoluition of biological tissue interface during tissue growth. Our model considers mechanical interactions between neighbouring cells and a secretion rate of new tissue material that is proportional to cell density. In this project we include cell proliferation, apoptosis (death) and embedment as stochastic processes. This model is solved using the `DifferentialEquations.jl` and `ElasticArrays.jl` packages. The solver uses `RK4()` or `Euler()` solving algorithms with a constant timestep $\Delta t$. The stochastic cell behaviour is implemented using a `PeriodicCallback` and occurs every $\delta t_{\text{event}}$ period of time. This package offers 1D and 2D simulations. In our case 1D simulations imply the evolution of a line segment which has periodic boundaries and 2D is a closed domain. In the case of 2D simulations you can choose an `inward` or `outward` growth simulation by setting the appropriate parameters.

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
given $f(\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert)$ is the restoration force function and $Vₙ$ is the velocity in the normal direction.


### Current models of Mechanical Relaxation and Normal Velocity
In this current version of the simulation code we use a nonlinear resorting force given by,

```math
f(\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert) = kₛ \bigg(\frac{1}{l_{0}} - \frac{1}{\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert }\bigg)
```
and a normal velocity which is proportional to density such that,
```math
Vₙ = k_{f}\rho_{i}, \hspace{0.5cm} \rho_{i} = \frac{1}{\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert}
```

