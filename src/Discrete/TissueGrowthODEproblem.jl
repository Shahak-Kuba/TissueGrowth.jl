using DifferentialEquations


"""
    ODE_fnc_1D!(du, u, p, t)

Define the ODE system for 1D mechanical relaxation with initial conditions. This function computes the derivatives `du` based on the current state `u` and parameters `p`. In this system there is a periodic boundary condition such that `x ∈ [0,2π]` with `i` spring boundary nodes and `u₋₁ = uₙ`, `uₙ₊₁ = u₁`

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

    uᵢ₊₁ = circshift(u',1)
    uᵢ₋₁ = circshift(u',-1)
    uᵢ₋₁[end,:] = uᵢ₋₁[end,:]+[dom;0]
    uᵢ₊₁[1,:] = uᵢ₊₁[1,:]-[dom;0]
    du .= ((1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁) +
                        Vₙ(uᵢ₋₁,u',uᵢ₊₁,kf,δt,growth_dir))'
    nothing
end


"""
    ODE_fnc_2D!(du, u, p, t)

Define the ODE system for 2D mechanical relaxation with initial conditions. This function computes the derivatives `du` based on the current state `u` and parameters `p`.

# Arguments
- `du`: Array to store the derivatives of `u`.
- `u`: Current position array of all spring boundary nodes.
- `p`: Parameters tuple `(N, kₛ, η, kf, l₀, δt, growth_dir)`.
- `t`: Current time.

# Description
Calculates the mechanical relaxation and normal velocity in a 1D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in `p`.
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
