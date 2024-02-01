"""
    PostCalcs1D(u, p)

Perform post-calculation for 1D simulation data.

This function computes various physical quantities like force, density, velocity, stress, and curvature for a given state `u` and parameters `p`.

# Arguments
- `u`: A state vector representing the positions of particles or cells.
- `p`: A tuple of parameters used in the calculations.

# Returns
A tuple containing the sum of forces, normal velocity, density, stress, and curvature for each element in the state vector.
"""
function PostCalcs1D(u, p)
    N, kₛ, η, kf, l₀, δt = p
    dom = 2*pi

    ∑F = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ = zeros(size(u, 1))
    ψ = zeros(size(u, 1))
    Κ = zeros(size(u, 1))

    uᵢ₊₁ = circshift(u',1)
    uᵢ₋₁ = circshift(u',-1)
    uᵢ₋₁[end,:] = uᵢ₋₁[end,:]+[dom;0]
    uᵢ₊₁[1,:] = uᵢ₊₁[1,:]-[dom;0]

    ∑F = diag((Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁)))
    density = ρ(uᵢ₊₁, u)
    ψ = ∑F ./ δ(uᵢ₊₁, u)
    Κ = κ(uᵢ₋₁,u,uᵢ₊₁)
    vₙx = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt,"2D")[:,1]
    vₙy = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt,"2D")[:,2]
    vₙ = .√(vₙx.^2 + vₙy.^2)

    return ∑F, vₙ, density, ψ, Κ
end


"""
    postSimulation1D(btype, sol, p)

Perform post simulation calculations for 1D simulation and return a comprehensive data structure with all relevant data.

This function processes the solution from a 1D simulation, calculating physical quantities and organizing them into a `SimResults_t` data structure.

# Arguments
- `btype`: The type of boundary condition or simulation.
- `sol`: The solution object from the simulation.
- `p`: Parameters used in the post calculations.

# Returns
An instance of `SimResults_t` containing the calculated data.
"""
function postSimulation1D(btype, sol, p)

    c = size(sol.t, 1)

    Area = Vector{Float64}(undef, c)
    Cell_Count = Vector{Float64}(undef, c)
    ∑F = Vector{Vector{Float64}}(undef, 0)
    ψ = Vector{Matrix{Float64}}(undef, 0)
    DENSITY = Vector{Matrix{Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Matrix{Float64}}(undef, 0)

    u = [(reshape(vec, 2, Int(length(vec)/2)))' for vec in sol.u]

    for ii in axes(u, 1)
        Area[ii] = Ω(u[ii]) # area calculation
        Cell_Count[ii] = size(u[ii],1)
        Fnet, nV, den, stre, kap = PostCalcs2D(u[ii], p)
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(DENSITY, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(btype, sol.t, u, ∑F, DENSITY, vₙ, Area, ψ, Κ, Cell_Count)
end


###################################################################################################

"""
    PostCalcs2D(u, p)

Perform post-calculation for 2D simulation data.

This function is similar to `PostCalcs1D`, but it is tailored for 2D simulation data. It calculates force, density, velocity, stress, and curvature in a 2D context.

# Arguments
- `u`: A 2D state vector representing the positions of particles or cells.
- `p`: A tuple of parameters used in the calculations.

# Returns
A tuple containing the sum of forces, normal velocity, density, stress, and curvature for each element in the state vector.
"""
function PostCalcs2D(u, p)
    N, kₛ, η, kf, l₀, δt = p

    #u = reshape(u, Int(length(u)/2), 2)

    ∑F = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ = zeros(size(u, 1))
    ψ = zeros(size(u, 1))
    Κ = zeros(size(u, 1))

    uᵢ₊₁ = circshift(u,1)
    uᵢ₋₁ = circshift(u,-1)

    ∑F = diag(Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) * transpose(τ(uᵢ₊₁,u))) + diag(Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) * transpose(τ(u,uᵢ₋₁)))
    #diag((Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁)))
    density = (ρ(uᵢ₊₁, u).+ρ(u, uᵢ₋₁))./2
    ψ = ∑F ./ δ(uᵢ₊₁, u)
    Κ = κ(uᵢ₋₁,u,uᵢ₊₁)
    vₙx = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt,"2D")[:,1]
    vₙy = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt,"2D")[:,2]
    vₙ = .√(vₙx.^2 + vₙy.^2)


    return ∑F, vₙ, density, ψ, Κ
end

"""
    postSimulation2D(btype, sol, p)

Perform post simulation calculations for 2D simulation and return a comprehensive data structure with all relevant data.

This function processes the solution from a 2D simulation, similarly to `postSimulation1D`, but adapted for 2D data.

# Arguments
- `btype`: The type of boundary condition or simulation.
- `sol`: The solution object from the simulation.
- `p`: Parameters used in the post calculations.

# Returns
An instance of `SimResults_t` containing the calculated data.
"""
function postSimulation2D(btype, sol, p)

    c = size(sol.t, 1)

    Area = Vector{Float64}(undef, c)
    Cell_Count = Vector{Float64}(undef, c)
    ∑F = Vector{Vector{Float64}}(undef, 0)
    ψ = Vector{Matrix{Float64}}(undef, 0)
    DENSITY = Vector{Matrix{Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Matrix{Float64}}(undef, 0)

    u = [(reshape(vec, 2, Int(length(vec)/2)))' for vec in sol.u]

    for ii in axes(u, 1)
        Area[ii] = Ω(u[ii]) # area calculation
        Cell_Count[ii] = size(u[ii],1)
        Fnet, nV, den, stre, kap = PostCalcs2D(u[ii], p)
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(DENSITY, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(btype, sol.t, u, ∑F, DENSITY, vₙ, Area, ψ, Κ, Cell_Count)
end

