# User defined force function

## Hookean Restoring force
hookean_restoring_force = (rᵢ, rⱼ, kₛ, l₀) -> kₛ .* ( δ(rⱼ,rᵢ) .- ones(size(rᵢ,1))*l₀ ) 
## Nonlinear restoring force
nonlinear_restoring_force = (rᵢ, rⱼ, kₛ, l₀) -> kₛ .* (ones(size(rᵢ,1),1) ./ l₀ - 1 ./ δ(rⱼ, rᵢ))
nonlinear_restoring_force2 = (rᵢ, rⱼ, kₛ, l₀) -> kₛ .* l₀.^2 .* (ones(size(rᵢ,1),1) ./ l₀ - 1 ./ δ(rⱼ, rᵢ))


# When changing force law make sure to run all of these
#FORCE_FNC = (rᵢ, rⱼ, kₛ, l₀) -> nonlinear_restoring_force(rᵢ, rⱼ, kₛ, l₀)

# Force functions used in ODEs
#Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) =  FORCE_FNC(rᵢ, rᵢ₊₁, kₛ, l₀) .* τ(rᵢ₊₁, rᵢ)
#Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = -FORCE_FNC(rᵢ, rᵢ₋₁, kₛ, l₀) .* τ(rᵢ, rᵢ₋₁)

function Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀, restoring_force)
    if restoring_force == "hookean"
        return hookean_restoring_force(rᵢ, rᵢ₊₁, kₛ, l₀).* τ(rᵢ₊₁, rᵢ)
    elseif restoring_force == "nonlinear"
        return nonlinear_restoring_force(rᵢ, rᵢ₊₁, kₛ, l₀).* τ(rᵢ₊₁, rᵢ)
    else
        return nonlinear_restoring_force2(rᵢ, rᵢ₊₁, kₛ, l₀).* τ(rᵢ₊₁, rᵢ)
    end
end

function Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀, restoring_force)
    if restoring_force == "hookean"
        return -hookean_restoring_force(rᵢ, rᵢ₋₁, kₛ, l₀).* τ(rᵢ, rᵢ₋₁)
    elseif restoring_force == "nonlinear"
        return -nonlinear_restoring_force(rᵢ, rᵢ₋₁, kₛ, l₀).* τ(rᵢ, rᵢ₋₁)
    else
        return -nonlinear_restoring_force2(rᵢ, rᵢ₋₁, kₛ, l₀).* τ(rᵢ, rᵢ₋₁)
    end
end

