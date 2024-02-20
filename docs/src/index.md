```@meta
CurrentModule = TissueGrowth
```

# TissueGrowth.jl

TissueGrowth.jl is a open source project package developed to simulate the evoluition of biological tissue interface during tissue growth. This package is used to solve both the system of ODEs defined in the discrete model and the PDE derived through taking the continuum limit of the discrete model.

# Discrete model description
Our model considers mechanical interactions between neighbouring cells and a secretion rate of new tissue material that is proportional to cell density. In this project we include cell proliferation, apoptosis (death) and embedment as stochastic processes. This model is solved using the `DifferentialEquations.jl` and `ElasticArrays.jl` packages. The solver uses `RK4()` or `Euler()` solving algorithms with a constant timestep $\Delta t$. The stochastic cell behaviour is implemented using a `PeriodicCallback` and occurs every $\delta t_{\text{event}}$ period of time. This package offers 1D and 2D simulations. In our case 1D simulations imply the evolution of a line segment which has periodic boundaries and 2D is a closed domain. In the case of 2D simulations you can choose an `inward` or `outward` growth simulation by setting the appropriate parameters.

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
f(\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert) = kₛ^{(m)} \bigg(\frac{1}{l_{0}} - \frac{1}{\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert }\bigg)
```
and a normal velocity which is proportional to density such that,
```math
Vₙ = k_{f}^{(m)}\rho_{i}, \hspace{0.5cm} \rho_{i} = \frac{1}{\lVert \mathbf{u}ᵢ₊₁ - \mathbf{u}ᵢ \rVert}
```

# Continuum model description
To derive the continuum limit of the discrete model we consider a very large number of springs per cell such that $m\rightarrow\infty$. There we are able to derive the evolution equation for spring density $\frac{\partial\rho}{\partial t}$ and then transform the result into an evolution equation for cell density $\frac{\partial q}{\partial t}$ using the relationship between spring and cell density $q(x,t) = \frac{\rho(x,t)}{m}$. Using appropriate parameter scalings we derive the evolution equation of cell density to be,
```math
\frac{\partial q}{\partial t} = -\frac{\partial^{2}\tilde{F}(q)}{\partial s^{2}} - qV(q)\kappa,
```
the mechanical relaxation is given by the first term and can be rewritten using the chain rule to look like a nonlinear diffusion,
```math
 -\frac{\partial^{2}\tilde{F}(q)}{\partial s^{2}} = \frac{\partial}{\partial s}\bigg(-\frac{\partial\tilde{F}(q)}{\partial q}\frac{\partial q}{\partial s}\bigg).
```
Given the choice of the nonlinear restoring force in the discrete models we retain a linear diffusion where,
```math
\tilde{F}(q) = D(q_{0} - q), \hspace{0.5cm} D = \frac{k_{s}^{(1)}}{\eta^{(1)}}q_{0}^{2},
```
```math
V(q) = k_{f}^{(1)}q,
```
$\frac{\partial^{2}}{\partial s^{2}}$ is the second derivative in terms of arc-length and $\bm{\gamma}(x,t)$ is an arbitrary parameterisation for the interface.

In cartesian form such that $\bm{\gamma} = (x,h(x))$ there is an extra term such that,

```math
\frac{\partial q}{\partial t} = -\frac{\partial^{2}\tilde{F}(q(x,t))}{\partial s^{2}} - q(x,t)V(q(x,t))\kappa - \bigg(\frac{\bm{\gamma}_{t}\cdot\bm{\hat{\tau}}}{\lVert \bm{\gamma}_{t} \rVert} \bigg)\frac{\partial q(x,t)}{\partial x}
```

For high to mid diffusivities i.e. $1 \geq D \geq 0.05$ we use a Finite Difference upwinding scheme. 

