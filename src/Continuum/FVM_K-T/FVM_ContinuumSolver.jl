include("FVM_SolverFncs.jl")

function FVM_SolveContinuumLim_Polar(D,kf,A,ρ₀,Tmax)
    # time and space discretisation
    T₀ = 0.0

    N = 10000
    Δt = (Tmax - T₀)/N
    t = LinRange(T₀,Tmax,Int64((Tmax - T₀)/Δt))
 
    m = 161
    Δθ = (2π)/m
    θ = Vector(LinRange(0.0,2π,m))
    pop!(θ)
    M = m - 1

    Φ = 1; # for minmod gradient 

    # allocating simulation memory
    ρ = zeros(N+1,M)
    r = zeros(N+1,M)
    dr = zeros(N+1,M)
    ddr = zeros(N+1,M)
    η = zeros(N+1,M)
    dη = zeros(N+1,M)
    σ = zeros(N+1,M)
    dσ = zeros(N+1,M)
    κ = zeros(N+1,M)

    rᵢ = zeros(size(r[1,:]))
    rᵢ₊₁ = zeros(size(r[1,:]))
    rᵢ₋₁ = zeros(size(r[1,:]))
    ηᵢ = zeros(size(η[1,:]))
    ηᵢ₊₁ = zeros(size(η[1,:]))
    ηᵢ₋₁ = zeros(size(η[1,:]))
    σᵢ = zeros(size(σ[1,:]))
    σᵢ₊₁ = zeros(size(σ[1,:]))
    σᵢ₋₁ = zeros(size(σ[1,:]))


    # setting up initial conditions
    # (interface: r₀)
    R = 3.0;
    for i in 1:M
        if θ[i] < π/2
            r[1,i] = min((R/cos(θ[i])), (R/sin(θ[i])))
        else
            r[1,i] = min((R/cos((θ[i]%(π/2)))), (R/sin((θ[i]%(π/2)))))
        end
    end
    # (σ₀)
    rᵢ .= r[1,:]
    rᵢ₊₁ .= circshift(rᵢ,1)
    rᵢ₋₁ .= circshift(rᵢ,-1)
    σ[1,:] .= (rᵢ₊₁.-rᵢ₋₁)./(2 .*Δt)

    # (η₀)
    η[1,:] .= ρ₀.*sqrt.(rᵢ.^2 + σᵢ.^2)
    
    # (ρ₀)
    ρ[1,:] .= ρ₀.*ones(size(ρ[1,:]))


    # simulation time loop 
    @time for n in 1:N

        rᵢ .= r[n,:]
        rᵢ₊₁ .= circshift(rᵢ,1)
        rᵢ₋₁ .= circshift(rᵢ,-1)
        σᵢ .= σ[n,:]
        σᵢ₊₁ .= circshift(σᵢ,1)
        σᵢ₋₁ .= circshift(σᵢ,-1)
        ηᵢ .= η[n,:]
        ηᵢ₊₁ .= circshift(ηᵢ,1)
        ηᵢ₋₁ .= circshift(ηᵢ,-1)

        # getting derivatives
        dr[n,:] .= minmod(rᵢ₋₁, rᵢ, rᵢ₊₁, Φ, Δθ)
        dσ[n,:] .= minmod(σᵢ₋₁, σᵢ, σᵢ₊₁, Φ, Δθ)
        dη[n,:] .= minmod(ηᵢ₋₁, ηᵢ, ηᵢ₊₁, Φ, Δθ)
        ddr[n,:] .= (rᵢ₋₁ .- 2 .*rᵢ .+ rᵢ₊₁)./(Δθ^2)

        # calculate curvature
        κ[n,:] .= -(rᵢ.^2 .- rᵢ.*ddr[n,:] .+ 2 .* dr[n,:])./((rᵢ.^2 .+ dr[n,:].^2).^1.5)

        # Finite Volume Method Kraganov-Tadmor scheme
        r⁺₁ = rᵢ .- (Δθ/2).*dr[n,:]
        r⁺₂ = rᵢ .+ (Δθ/2).*dr[n,:]
        σ⁺₁ = σᵢ .- (Δθ/2).*dσ[n,:]
        σ⁺₂ = σᵢ .+ (Δθ/2).*dσ[n,:]
        η⁺₁ = ηᵢ .- (Δθ/2).*dη[n,:]
        η⁺₂ = ηᵢ .+ (Δθ/2).*dη[n,:]

        r⁻₁ = rᵢ .- (Δθ/2).*dr[n,:]
        r⁻₂ = rᵢ .+ (Δθ/2).*dr[n,:]
        σ⁻₁ = σᵢ .- (Δθ/2).*dσ[n,:]
        σ⁻₂ = σᵢ .+ (Δθ/2).*dσ[n,:]
        η⁻₁ = ηᵢ .- (Δθ/2).*dη[n,:]
        η⁻₂ = ηᵢ .+ (Δθ/2).*dη[n,:]

        # derivatives for Kraganov-Tadmor scheme
        dσ⁺₁ = minmod(circshift(σ⁺₁,-1), σ⁺₁, circshift(σ⁺₁,1), Φ, Δθ)
        dσ⁺₂ = minmod(circshift(σ⁺₂,-1), σ⁺₂, circshift(σ⁺₂,1), Φ, Δθ)
        dη⁺₁ = minmod(circshift(η⁺₁,-1), η⁺₁, circshift(η⁺₁,1), Φ, Δθ)
        dη⁺₂ = minmod(circshift(η⁺₂,-1), η⁺₂, circshift(η⁺₂,1), Φ, Δθ)

        dσ⁻₁ = minmod(circshift(σ⁻₁,-1), σ⁻₁, circshift(σ⁻₁,1), Φ, Δθ)
        dσ⁻₂ = minmod(circshift(σ⁻₂,-1), σ⁻₂, circshift(σ⁻₂,1), Φ, Δθ)
        dη⁻₁ = minmod(circshift(η⁻₁,-1), η⁻₁, circshift(η⁻₁,1), Φ, Δθ)
        dη⁻₂ = minmod(circshift(η⁻₂,-1), η⁻₂, circshift(η⁻₂,1), Φ, Δθ)

        β⁺₁ = β(r⁺₁, σ⁺₁, η⁺₁, kf, D, dη⁺₁, dσ⁺₁)
        γ⁺₁ = γ(r⁺₁, σ⁺₁, η⁺₁, kf, D, dσ⁺₁)
        β⁺₂ = β(r⁺₂, σ⁺₂, η⁺₂, kf, D, dη⁺₂, dσ⁺₂)
        γ⁺₂ = γ(r⁺₂, σ⁺₂, η⁺₂, kf, D, dσ⁺₂)

        λ⁺₁ = λ( β⁺₁, γ⁺₁, r⁺₁)
        λ⁺₂ = λ( β⁺₂, γ⁺₂, r⁺₂)
        a⁺val = a⁺(λ⁺₁, λ⁺₂)

        β⁻₁ = β(r⁻₁, σ⁻₁, η⁻₁, kf, D, dη⁻₁, dσ⁻₁)
        γ⁻₁ = γ(r⁻₁, σ⁻₁, η⁻₁, kf, D, dσ⁻₁)
        β⁻₂ = β(r⁻₂, σ⁻₂, η⁻₂, kf, D, dη⁻₂, dσ⁻₂)
        γ⁻₂ = γ(r⁻₂, σ⁻₂, η⁻₂, kf, D, dσ⁻₂)

        λ⁻₁ = λ( β⁻₁, γ⁻₁, r⁻₁)
        λ⁻₂ = λ( β⁻₂, γ⁻₂, r⁻₂)
        a⁻val = a⁻(λ⁻₁, λ⁻₂)
        
        κ⁺₁ = -(r⁺₁.^2 .- r⁺₁.*dσ⁺₁ .+ 2 .* σ⁺₁.^2)./((r⁺₁.^2 .+ σ⁺₁.^2).^1.5)
        κ⁺₂ = -(r⁺₂.^2 .- r⁺₂.*dσ⁺₂ .+ 2 .* σ⁺₂.^2)./((r⁺₂.^2 .+ σ⁺₂.^2).^1.5)
        κ⁻₁ = -(r⁻₁.^2 .- r⁻₁.*dσ⁻₁ .+ 2 .* σ⁻₁.^2)./((r⁻₁.^2 .+ σ⁻₁.^2).^1.5)
        κ⁻₂ = -(r⁻₂.^2 .- r⁻₂.*dσ⁻₂ .+ 2 .* σ⁻₂.^2)./((r⁻₂.^2 .+ σ⁻₂.^2).^1.5)

        F₁⁺₁ = f₁(r⁺₁,σ⁺₁,η⁺₁)
        F₁⁺₂ = f₁(r⁺₂,σ⁺₂,η⁺₂)
        F₁⁻₁ = f₁(r⁻₁,σ⁻₁,η⁻₁)
        F₁⁻₂ = f₁(r⁻₂,σ⁻₂,η⁻₂)

        F₂⁺₁ = f₂(r⁺₁,σ⁺₁,η⁺₁,kf)
        F₂⁺₂ = f₂(r⁺₂,σ⁺₂,η⁺₂,kf)
        F₂⁻₁ = f₂(r⁻₁,σ⁻₁,η⁻₁,kf)
        F₂⁻₂ = f₂(r⁻₂,σ⁻₂,η⁻₂,kf) 

        F₃⁺₁ = f₃(r⁺₁,σ⁺₁,η⁺₁,kf,D,dη⁺₁,κ⁺₁)
        F₃⁺₂ = f₃(r⁺₂,σ⁺₂,η⁺₂,kf,D,dη⁺₂,κ⁺₂)
        F₃⁻₁ = f₃(r⁻₁,σ⁻₁,η⁻₁,kf,D,dη⁻₁,κ⁻₁)
        F₃⁻₂ = f₃(r⁻₂,σ⁻₂,η⁻₂,kf,D,dη⁻₂,κ⁻₂)

        H⁺₁ = 0.5.*(F₁⁺₁ .+ F₁⁺₂) .- ((a⁺val./2).*(r⁺₁ .- r⁺₂))
        H⁻₁ = 0.5.*(F₁⁻₁ .+ F₁⁻₂) .- ((a⁻val./2).*(r⁻₁ .- r⁻₂))

        H⁺₂ = 0.5.*(F₂⁺₁ .+ F₂⁺₂) .- ((a⁺val./2).*(σ⁺₁ .- σ⁺₂))
        H⁻₂ = 0.5.*(F₂⁻₁ .+ F₂⁻₂) .- ((a⁻val./2).*(σ⁻₁ .- σ⁻₂))

        H⁺₃ = 0.5.*(F₃⁺₁ .+ F₃⁺₂) .- ((a⁺val./2).*(η⁺₁ .- η⁺₂))
        H⁻₃ = 0.5.*(F₃⁻₁ .+ F₃⁻₂) .- ((a⁻val./2).*(η⁻₁ .- η⁻₂))

        L₁ = -(1/Δθ).*(H⁺₁ .- H⁻₁) .- kf.*(ηᵢ./rᵢ)
        L₂ = -(1/Δθ).*(H⁺₂ .- H⁻₂)
        L₃ = -(1/Δθ).*(H⁺₃ .- H⁻₃) .- A.*ηᵢ

        r[n+1,:] .= r[n,:] + Δt.*L₁
        σ[n+1,:] .= σ[n,:] + Δt.*L₂
        η[n+1,:] .= η[n,:] + Δt.*L₃
        ρ[n+1,:] .= η[n+1,:]./sqrt.(r[n+1,:].^2 + σ[n+1,:].^2)

    end
    return r,σ,η,ρ,θ
end 