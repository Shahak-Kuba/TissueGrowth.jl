function PostCalcs1D(u, p)
    N, kₛ, η, kf, l₀, δt = p
    dom = 2*pi

    ∑F = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ = zeros(size(u, 1))
    ψ = zeros(size(u, 1))
    Κ = zeros(size(u, 1))

    uᵢ₊₁ = circshift(u,-1)
    uᵢ₋₁ = circshift(u,1)
    uᵢ₋₁[1,:] = uᵢ₋₁[1,:]-[dom;0]
    uᵢ₊₁[end,:] = uᵢ₊₁[end,:]+[dom;0]
    ∑F = diag((Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀)) * transpose(τ(uᵢ₊₁,uᵢ₋₁)))
    density = ρ(uᵢ₊₁, u)
    ψ = ∑F ./ δ(uᵢ₊₁, u)
    Κ = κ(uᵢ₋₁,u,uᵢ₊₁)
    vₙx = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt,"1D")[:,1]
    vₙy = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kf,δt,"1D")[:,2]
    vₙ = .√(vₙx.^2 + vₙy.^2)

    return ∑F, vₙ, density, ψ, Κ
end


""" 
postSimulation1D()

Function to perform post simulation calculations which returns a data structure which contains all data
"""

function postSimulation1D(btype, sol, p)

    c = size(sol.t, 1)

    Area = Vector{Float64}(undef, c)
    ∑F = Vector{Vector{Float64, Vector{Float64}}}(undef, 0)
    ψ = Vector{Matrix{Float64, Vector{Float64}}}(undef, 0)
    DENSITY = Vector{ElasticMatrix{Float64, Vector{Float64}}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Matrix{Float64, Vector{Float64}}}(undef, 0)

    
    for ii in axes(sol.u, 1)
        Area[ii] = Ω(sol.u[ii]) # area calculation
        #append!(sol.u[ii], sol.u[ii][:,1]) # closing the domain Ω
        Fnet, nV, den, stre, kap = PostCalcs1D(sol.u[ii], p)
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(DENSITY, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(btype, sol.t, sol.u, ∑F, DENSITY, vₙ, Area, ψ, Κ)
end


###################################################################################################


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
postSimulation2D()

Function to perform post simulation calculations which returns a data structure which contains all data
"""

function postSimulation2D(btype, sol, p)

    c = size(sol.t, 1)

    Area = Vector{Float64}(undef, c)
    ∑F = Vector{Vector{Float64}}(undef, 0)
    ψ = Vector{Matrix{Float64}}(undef, 0)
    DENSITY = Vector{Matrix{Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Matrix{Float64}}(undef, 0)

    u = [(reshape(vec, 2, Int(length(vec)/2)))' for vec in sol.u]

    for ii in axes(u, 1)
        Area[ii] = Ω(u[ii]) # area calculation
        #append!(sol.u[ii], sol.u[ii][:,1]) # closing the domain Ω
        Fnet, nV, den, stre, kap = PostCalcs2D(u[ii], p)
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(DENSITY, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(btype, sol.t, u, ∑F, DENSITY, vₙ, Area, ψ, Κ)
end

function postSimulation2D_Elastic(btype, sol, p)

    c = size(sol.t, 1)

    Area = Vector{Float64}(undef, c)
    ∑F = Vector{ElasticVector{Float64,Vector{Float64}}}(undef, 0)
    ψ = Vector{ElasticMatrix{Float64,Vector{Float64}}}(undef, 0)
    DENSITY = Vector{ElasticMatrix{Float64,Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{ElasticMatrix{Float64,Vector{Float64}}}(undef, 0)


    for ii in axes(u, 1)
        Area[ii] = Ω(u[ii]) # area calculation
        #append!(sol.u[ii], sol.u[ii][:,1]) # closing the domain Ω
        Fnet, nV, den, stre, kap = PostCalcs2D(u[ii], p)
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(DENSITY, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_Elastic_t(btype, sol.t, u, ∑F, DENSITY, vₙ, Area, ψ, Κ)
end