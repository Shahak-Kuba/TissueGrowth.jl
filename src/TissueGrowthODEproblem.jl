using DifferentialEquations

# Mechanical Relaxation ODE problem
# Open boundary ODE Function (PERIODIC Boundary)
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


function ODE_fnc_2D_init!(du,u,p,t) 
    N,kₛ,η,kf,l₀,δt,growth_dir = p
    uᵢ₊₁ = circshift(u',1)
    uᵢ₋₁ = circshift(u',-1)
    du .= (1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁)
    nothing
end

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
SetupODEproblem(btype,N,R₀,kₛ,η,kf,l₀,δt,Tmax) is a function to create the split ODE problem, for DifferentalEquations.jl. 

The inputs to this equations include:

btype: boundary type, 

N: number of springs in the system, 

R₀: radius for a circle based problem (this will get adjusted for different shapes), 

kₛ,η: the prescaled mechanical relaxation coefficients, 

kf: the amount of tissue produced per cell per unit time, 

l₀: resting spring length, 

δt: Euler timestep size.

Tmax: end of simulation time.

growth_dir: direction of tissue growth. Options: "inward", "outward". For 1D simulations using outward is growth in positive y direction.
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