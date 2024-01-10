using LinearAlgebra

"""
    δ(rᵢ₊₁, rᵢ)

Calculate the Euclidean distance between two points `rᵢ₊₁` and `rᵢ`.

# Arguments
- `rᵢ₊₁`: The first point in space.
- `rᵢ`: The second point in space.

# Returns
The Euclidean distance between the two points.
"""
δ(rᵢ₊₁, rᵢ) = .√(sum((rᵢ₊₁ - rᵢ).^2,dims=size(rᵢ)))

"""
    τ(rᵢ₊₁, rᵢ₋₁)

Calculate the unit tangent vector between two neighboring points `rᵢ₊₁` and `rᵢ₋₁`.

# Arguments
- `rᵢ₊₁`: The point after the central point in space.
- `rᵢ₋₁`: The point before the central point in space.

# Returns
The unit tangent vector between the two points.
"""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ - rᵢ₋₁) ./ δ(rᵢ₊₁, rᵢ₋₁)

"""
    n(rᵢ₊₁, rᵢ₋₁, type)

Calculate the unit normal vector at a point `rᵢ` between two neighboring points `rᵢ₊₁` and `rᵢ₋₁`. The orientation of the normal vector depends on the specified `type`.

# Arguments
- `rᵢ₊₁`: The point after the central point in space.
- `rᵢ₋₁`: The point before the central point in space.
- `type`: A string specifying the orientation of the normal vector, either "inward" or any other value for outward orientation.

# Returns
The unit normal vector at the point.
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
    Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)

Calculate the spring force (Nonlinear) for mechanical relaxation in the positive direction.

# Arguments
- `rᵢ`: The current point in space.
- `rᵢ₊₁`: The point after the current point in space.
- `rᵢ₋₁`: The point before the current point in space.
- `kₛ`: Spring coefficient.
- `l₀`: Resting length of the spring.

# Returns
The spring force in the positive direction.
"""
Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = kₛ .* l₀.^2 .* (ones(size(rᵢ,1),1) ./ l₀ - 1 ./ δ(rᵢ₊₁, rᵢ)) .* τ(rᵢ₊₁, rᵢ)



"""
    Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀)

Calculate the spring force (Nonlinear) for mechanical relaxation in the negative direction.

# Arguments
- `rᵢ`: The current point in space.
- `rᵢ₊₁`: The point after the current point in space.
- `rᵢ₋₁`: The point before the current point in space.
- `kₛ`: Spring coefficient.
- `l₀`: Resting length of the spring.

# Returns
The spring force in the negative direction.
"""
Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = -kₛ .* l₀.^2 .* (ones(size(rᵢ,1),1) ./ l₀ - 1 ./ δ(rᵢ, rᵢ₋₁)) .* τ(rᵢ, rᵢ₋₁)


"""
    ρ(rᵢ₊₁, rᵢ)

Calculate the reciprocal of the distance (interpreted as density) between two points `rᵢ₊₁` and `rᵢ`.

# Arguments
- `rᵢ₊₁`: The first point in space.
- `rᵢ`: The second point in space.

# Returns
The reciprocal of the distance between the two points.
"""
ρ(rᵢ₊₁, rᵢ) = 1 ./ δ(rᵢ₊₁, rᵢ);


"""
    Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)

Calculate the normal velocity of the interface such that `Vₙ` is proportional to `ρ`. The direction of the normal vector is determined by `type`.

# Arguments
- `rᵢ₋₁`: The point before the current point in space.
- `rᵢ`: The current point in space.
- `rᵢ₊₁`: The point after the current point in space.
- `kf`: The amount of tissue produced per unit area per unit time.
- `δt`: The time step.
- `type`: A string specifying the orientation of the normal vector, either "inward" or any other value for outward orientation.

# Returns
The normal velocity of the interface.
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
    κ(rᵢ₋₁, rᵢ, rᵢ₊₁)

Approximate the curvature of a shape using the Menger method. This method is based on the areas of triangles formed by consecutive triplets of points.

# Arguments
- `rᵢ₋₁`: The point before the current point in space.
- `rᵢ`: The current point in space.
- `rᵢ₊₁`: The point after the current point in space.

# Returns
The approximated curvature at the point `rᵢ`.

# Reference
Anoshkina, Elena V., Alexander G. Belyaev, and Hans-Peter Seidel. "Asymtotic Analysis of Three-Point Approximations of Vertex Normals and Curvatures." VMV. 2002.
"""
function κ(rᵢ₋₁, rᵢ, rᵢ₊₁)

    A = ωκ(rᵢ₋₁,rᵢ,rᵢ₊₁)
    l1 = δ(rᵢ₋₁,rᵢ)
    l2 = δ(rᵢ,rᵢ₊₁)
    l3 = δ(rᵢ₋₁,rᵢ₊₁)

    return (4*A)./(l1.*l2.*l3)
end


