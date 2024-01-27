
function minmod(varᵢ₋₁, varᵢ, varᵢ₊₁, Φ, Δx)
    # flux limiter function given by [Kraganov-Tadmor 2000] FVM Scheme equation (4.9), 
    # this function is a vectorised code based of the implementation seen in Alias & Buenzli 2017.
    s1 = sign.(Φ.*(varᵢ₊₁ .- varᵢ)./Δx)
    s2 = sign.((varᵢ₊₁ .- varᵢ₋₁)/(2 .*Δx))
    s3 = sign.(Φ.*(varᵢ .- varᵢ₋₁)./Δx)

    output = zeros(size(varᵢ))

    for i in eachindex(output)
        if s1[i] == s2[i] && s1[i] == s3[i]
            output[i] = s1[i]*min(abs(Φ*(varᵢ₊₁[i] - varᵢ[i])/Δx), abs((varᵢ₊₁[i] - varᵢ₋₁[i])/(2 *Δx)), abs(Φ*(varᵢ[i] - varᵢ₋₁[i])/Δx))
        end
    end

    return output
end


# Functions from Modified conservation laws
function f₁(r, σ, η)
    return zeros(size(r))
end

function f₂(r, σ, η, kf)
    return kf.*η./r
end

function f₃(r, σ, η, kf, D, dη, κ)
    return (kf.*σ.*(η.^2))./(r.^3 + r.*σ.^2) .- (D.*dη)./(r.^2 .+ σ.^2) .+ (D.*η.*κ.*σ)./(r.*sqrt.(r.^2 .+ σ.^2)) .+ (2 .*D.*η.*σ)./(r.^3 + r.*σ.^2)
end


# Functions to calculate eigenvalues of the Jacobian
function β(r, σ, η, kf, D, dη, dσ)
    return (kf.*(η.^2).*r)./((r.^2 .+ σ.^2).^2) .- (kf.*(η.^2).*(σ.^2))./(r.*(r.^2 .+ σ.^2).^2) .+
            (D.*η.*r)./((r.^2 .+ σ.^2).^2) .- (8 .*D.*η.*(σ.^2))./(r.*(r.^2 .+ σ.^2).^2) .+
            (2 .*D.*σ.*dη)./((r.^2 .+ σ.^2).^2) .+ (4 .*D.*η.*(σ.^2).*r)./((r.^2 .+ σ.^2).^3) .+
            (8 .*D.*η.*(σ.^4))./(r.*(r.^2 .+ σ.^2).^3) .+ (D.*η.*dσ)./((r.^2 .+ σ.^2).^2) .-
            (4 .*D.*η.*(σ.^2).*dσ)./((r.^2 .+ σ.^2).^3)
end

function γ(r, σ, η, kf, D, dσ)
    return (2 .*kf.*σ.*η .+ 2 .*D.*σ)./(r.*(r.^2 .+ σ.^2)) .- 
           (D.*σ.*r)./((r.^2 .+ σ.^2).^2) .- 
           (2 .*D.*(σ.^3))./(r.*(r.^2 .+ σ.^2).^2) .+
           (D.*σ.*dσ)./((r.^2 .+ σ.^2).^2)
end

# Solving and returning max eigenvalues (i.e. roots to [-1 γ β/r 0] and characteristic equation given by "-λ³ + γλ² + (β/r)λ = 0" )
function λ(Β,Γ,r)
    λ_max = zeros(size(r));
    λ₊ = real.(abs.((Γ .+ sqrt.(Complex.(Γ.^2 .- 4 .*(Β./r))))./2))
    λ₋ = real.(abs.((Γ .- sqrt.(Complex.(Γ.^2 .- 4 .*(Β./r))))./2))
    for i in eachindex(λ_max)
        λ_max[i] = max(λ₊[i], λ₋[i], 0.0)
    end
    return λ_max
end

function a⁺(λ⁺₁, λ⁺₂)
    return max.(λ⁺₁, λ⁺₂)
end

function a⁻(λ⁻₁, λ⁻₂)
    return max.(λ⁻₁, λ⁻₂)
end