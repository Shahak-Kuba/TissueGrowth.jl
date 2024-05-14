
"""
    SetupODEproblem(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, prolif, death, embed, α, β, γ, dist_type)

Set up and configure a 1D or 2D ODE problem for mechanical relaxation simulations in tissue growth.

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
function SetupODEproblem(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,growth_dir,domain_type,prolif,death,embed,β,γ,Ot,dist_type,restoring_force)
    l₀ = l₀/m
    kₛ = kₛ*m
    η = η/m
    kf = kf/m
    u0 = u0SetUp(btype,R₀,M,dist_type,domain_type)
    p = (m,kₛ,η,kf,l₀,δt,growth_dir,domain_type,btype,restoring_force,prolif,death,embed,β,γ,Ot)
    tspan = (0.0,Tmax)
    return ODEProblem(Growth_ODE2!,u0,tspan,p), p
end

# Setup for when cell density limit is applied
function SetupODEproblem(btype,M,m,R₀,kₛ,η,kf,l₀,δt,Tmax,growth_dir,domain_type,prolif,death,embed,β,γ,Ot,dist_type,restoring_force,q_lim)
    l₀ = l₀/m
    kₛ = kₛ*m
    η = η/m
    kf = kf/m
    u0 = u0SetUp(btype,R₀,M,dist_type,domain_type)
    p = (m,kₛ,η,kf,l₀,δt,growth_dir,domain_type,btype,restoring_force,prolif,death,embed,β,γ,Ot, q_lim)
    tspan = (0.0,Tmax)
    return ODEProblem(Growth_ODE2!,u0,tspan,p), p
end