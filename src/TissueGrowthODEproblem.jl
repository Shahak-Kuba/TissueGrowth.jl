using DifferentialEquations

"""
    ODE_fnc_1D_init!(du, u, p, t)

Define the ODE system for 1D mechanical relaxation with initial conditions. This function computes the derivatives `du` based on the current state `u` and parameters `p`.

# Arguments
- `du`: Array to store the derivatives of `u`.
- `u`: Current state array.
- `p`: Parameters tuple `(N, kₛ, η, kf, l₀, δt, growth_dir)`.
- `t`: Current time.

# Description
Calculates the mechanical relaxation in a 1D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in `p`. This ODE problem can be used to mechanically relax the system to have equal cell densities when including normal velocity in simulation.
"""

function ODE_fnc_1D_init!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt,growth_dir = p
    dom = 2*pi;
    uᵢ₊₁ = circshift(u,-1)
    uᵢ₋₁ = circshift(u,1)
    uᵢ₋₁[:,1] = uᵢ₋₁[:,1]-[dom;0]
    uᵢ₊₁[:,end] = uᵢ₊₁[:,end]+[dom;0]
    du .= (1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁)
    nothing
end

"""
    ODE_fnc_1D!(du, u, p, t)

Define the ODE system for 1D mechanical relaxation with initial conditions. This function computes the derivatives `du` based on the current state `u` and parameters `p`.

# Arguments
- `du`: Array to store the derivatives of `u`.
- `u`: Current state array.
- `p`: Parameters tuple `(N, kₛ, η, kf, l₀, δt, growth_dir)`.
- `t`: Current time.

# Description
Calculates the mechanical relaxation and normal velocity in a 1D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in `p`.
"""
function ODE_fnc_1D!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt,growth_dir = p
    dom = 2*pi;
    uᵢ₊₁ = circshift(u',-1)
    uᵢ₋₁ = circshift(u',1)
    uᵢ₋₁[:,1] = uᵢ₋₁[:,1]-[dom;0]
    uᵢ₊₁[:,end] = uᵢ₊₁[:,end]+[dom;0]
    du .= ((1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁) +
                       Vₙ(uᵢ₋₁,u',uᵢ₊₁,kf,δt,growth_dir))'
   
    nothing
end

"""
    ODE_fnc_2D_init!(du, u, p, t)

Define the ODE system for 1D mechanical relaxation with initial conditions. This function computes the derivatives `du` based on the current state `u` and parameters `p`.

# Arguments
- `du`: Array to store the derivatives of `u`.
- `u`: Current state array.
- `p`: Parameters tuple `(N, kₛ, η, kf, l₀, δt, growth_dir)`.
- `t`: Current time.

# Description
Calculates the mechanical relaxation in a 2D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in `p`. This ODE problem can be used to mechanically relax the system to have equal cell densities when including normal velocity in simulation.
"""
function ODE_fnc_2D_init!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt,growth_dir = p
    uᵢ₊₁ = circshift(u',1)
    uᵢ₋₁ = circshift(u',-1)
    du .= (1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁)
    nothing
end

"""
    ODE_fnc_1D!(du, u, p, t)

Define the ODE system for 2D mechanical relaxation with initial conditions. This function computes the derivatives `du` based on the current state `u` and parameters `p`.

# Arguments
- `du`: Array to store the derivatives of `u`.
- `u`: Current state array.
- `p`: Parameters tuple `(N, kₛ, η, kf, l₀, δt, growth_dir)`.
- `t`: Current time.

# Description
Calculates the mechanical relaxation and normal velocity in a 1D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in `p`. Using the discrete equation:
```math
duᵢ/dt = 1/η(Fₛ⁺ - Fₛ⁻)τ + Vₙn 
```
"""
function ODE_fnc_2D!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt,growth_dir = p
    #u_mat = reshape(u, 2, Int(length(u)/2))'
    uᵢ₊₁ = circshift(u',1)
    uᵢ₋₁ = circshift(u',-1)
    du .= ((1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁) +
                       Vₙ(uᵢ₋₁,u',uᵢ₊₁,kf,δt,growth_dir))'
    nothing
end


"""
    SetupODEproblem1D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, dist_type)

Create and configure a 1D ODE problem for mechanical relaxation simulations.

This function sets up the initial conditions and parameters for a 1D ODE problem based on the specified boundary type and other physical parameters.

# Arguments
- `btype`: Boundary type (e.g., 'circle', 'triangle').
- `M`, `m`: Number of springs/particles in the system.
- `R₀`: Initial radius for circular boundary problems.
- `kₛ`, `η`: Prescaled mechanical relaxation coefficients.
- `kf`: Tissue production rate per cell per unit time.
- `l₀`: Resting spring length.
- `δt`: Euler timestep size.
- `Tmax`: End of simulation time.
- `growth_dir`: Direction of tissue growth ('inward' or 'outward').
- `dist_type`: Distribution type for node placement.

# Returns
Configured ODEProblem instance and parameters for the 1D simulation.
"""
function SetupODEproblem1D(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,growth_dir,dist_type)
    l₀ = l₀/m
    kₛ = kₛ*m
    kf = kf/m
    η = η/m
    # setting up initial conditions
    u0 = u0SetUp(btype,R₀,M,dist_type)
    # solving ODE problem
    p = (M,kₛ,η,kf,l₀,δt,growth_dir)
    tspan = (0.0,Tmax)
    return ODEProblem(ODE_fnc_1D!,u0,tspan,p), p
end

"""
    SetupODEproblem2D(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, prolif, death, embed, α, β, γ, dist_type)

Set up and configure a 2D ODE problem for mechanical relaxation simulations in tissue growth.

This function initializes the conditions and parameters for a 2D ODE problem based on the specified boundary type, physical parameters, and cell behaviors. It then constructs an ODEProblem object, ready for solving with DifferentialEquations.jl.

# Arguments
- `btype`: Boundary type (e.g., 'circle', 'triangle').
- `M`: Total number of springs along the interface.
- `m`: Number of springs per cell.
- `R₀`: Initial radius or characteristic length of the shape.
- `kₛ`: Spring stiffness coefficient.
- `η`: Viscosity or damping coefficient.
- `kf`: Tissue production rate per cell.
- `l₀`: Resting length of the spring per cell.
- `δt`: Time step for the numerical integration.
- `Tmax`: Total simulation time.
- `growth_dir`: Direction of tissue growth ('inward' or 'outward').
- `prolif`: Boolean flag for cell proliferation.
- `death`: Boolean flag for cell death.
- `embed`: Boolean flag for cell embedding.
- `α, β, γ`: Parameters for the cell behaviors.
- `dist_type`: Distribution type for node placement.

# Returns
- `ODEProblem`: An ODE problem instance set up with the specified parameters and initial conditions.
- `p`: A tuple containing the parameters used in setting up the ODE problem.
"""
function SetupODEproblem2D(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,growth_dir,prolif,death,embed,α,β,γ,dist_type)
    l₀ = l₀/m
    kₛ = kₛ*m
    η = η/m
    kf = kf/m
    u0 = u0SetUp(btype,R₀,M,dist_type)
    #plotInitialCondition(u0)
    # solving ODE problem
    p = (m,kₛ,η,kf,l₀,δt,growth_dir,prolif,death,embed,α,β,γ)
    tspan = (0.0,Tmax)
    return ODEProblem(ODE_fnc_2D!,u0,tspan,p), p
end