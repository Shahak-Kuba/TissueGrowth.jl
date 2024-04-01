## Nonlinear restoring force

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
Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = -(kₛ .* l₀.^2 .* (ones(size(rᵢ,1),1) ./ l₀ - 1 ./ δ(rᵢ, rᵢ₋₁)) .* τ(rᵢ, rᵢ₋₁))