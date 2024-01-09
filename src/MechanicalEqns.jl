using LinearAlgebra

"""
δ is the length between two nodes in the positive direction
"""
δ(rᵢ₊₁, rᵢ) = .√(sum((rᵢ₊₁ - rᵢ).^2,dims=size(rᵢ)))

"""
τ calculates the unit tangent vector between two neighbouring points rᵢ₊₁ and rᵢ₋₁ of a central point rᵢ
"""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ - rᵢ₋₁) ./ δ(rᵢ₊₁, rᵢ₋₁)

"""
n calculates the unit normal vector between two neighbouring points rᵢ₊₁ and rᵢ₋₁ of a central point rᵢ
"""
function n(rᵢ₊₁, rᵢ₋₁,type) 
    if type == "inward"
        -oftype(τ(rᵢ₊₁, rᵢ₋₁),vcat(transpose.([-τ(rᵢ₊₁, rᵢ₋₁)[:,2], τ(rᵢ₊₁, rᵢ₋₁)[:,1]])...)')
    else
        oftype(τ(rᵢ₊₁, rᵢ₋₁),vcat(transpose.([-τ(rᵢ₊₁, rᵢ₋₁)[:,2], τ(rᵢ₊₁, rᵢ₋₁)[:,1]])...)')
    end
end
        #n(rᵢ₊₁, rᵢ₋₁) = (τv = τ(rᵢ₊₁, rᵢ₋₁);return (-τv[2], τv[1]))
"""
Fₛ⁺ is the spring force (Nonlinear) used for mechanical relaxation in the positive direction.
l = length of the spring
l₀ = resting length of the spring
kₛ = spring coefficient
"""

Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = kₛ .* l₀.^2 .* (ones(size(rᵢ,1),1) ./ l₀ - 1 ./ δ(rᵢ₊₁, rᵢ)) .* τ(rᵢ₊₁, rᵢ)

"""
Fₛ⁻ is the spring force (Nonlinear) used for mechanical relaxation in the negative direction.
l = length of the spring
l₀ = resting length of the spring
kₛ = spring coefficient
"""

Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = -kₛ .* l₀.^2 .* (ones(size(rᵢ,1),1) ./ l₀ - 1 ./ δ(rᵢ, rᵢ₋₁)) .* τ(rᵢ, rᵢ₋₁)

"""

"""
ρ(rᵢ₊₁, rᵢ) = 1 ./ δ(rᵢ₊₁, rᵢ);

"""
Vₙ is the normal velocity of the interface such that Vₙ ∝ ρ
ρ = the density of the cell
kf = the amount of tissue produced per unit area per unit time
"""

function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt,type)
    ρₗ = ρ(rᵢ, rᵢ₋₁)
    ρᵣ = ρ(rᵢ₊₁, rᵢ)
    Vₗ = kf .* ρₗ
    Vᵣ = kf .* ρᵣ

    nₗ = n(rᵢ₋₁, rᵢ,type)
    nᵣ = n(rᵢ, rᵢ₊₁,type)

    rₘ₁ = rᵢ - (Vₗ .* nₗ .* δt)
    rₗ = rᵢ₋₁ - (Vₗ .* nₗ .* δt)
    rₘ₂ = rᵢ - (Vᵣ .* nᵣ .* δt)
    rᵣ = rᵢ₊₁ - (Vᵣ .* nᵣ .* δt)

    return (lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ) - rᵢ) ./ δt
end


"""
κ(rᵢ₋₁,rᵢ,rᵢ₊₁) approximated the curvature of the shape using Menger method to approximate curvature
"""
# Anoshkina, Elena V., Alexander G. Belyaev, and Hans-Peter Seidel. "Asymtotic Analysis of Three-Point Approximations of Vertex Normals and Curvatures." VMV. 2002.

function κ(rᵢ₋₁, rᵢ, rᵢ₊₁)

    A = ωκ(rᵢ₋₁,rᵢ,rᵢ₊₁)
    l1 = δ(rᵢ₋₁,rᵢ)
    l2 = δ(rᵢ,rᵢ₊₁)
    l3 = δ(rᵢ₋₁,rᵢ₊₁)

    return (4*A)./(l1.*l2.*l3)
end


